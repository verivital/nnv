%% Robustness verification of a NN (L infinity adversarial attack)
%  if f(x) = y, then forall x' in X s.t. ||x - x'||_{\infty} <= eps,
%  then f(x') = y = f(x)
%
% certified robustness: 
%  what % of correctly classified inputs satisfies the above provably

%% Load data into NNV

% Load network 
mnist_model = load('mnist_model.mat');

% Create NNV model
net = matlab2nnv(mnist_model.net);

% Load data (no download necessary)
digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
    'nndatasets','DigitDataset');
% Images
imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

% Shuffle images for evaluation
rng(0); % set random seed for reproducibility
imds = shuffle(imds);

% Get a smaller subset for timing constraints
nClass = 100; % get 100 image per class (10% of whole dataset)
[imdsEval, ~] = splitEachLabel(imds, 100);

N = length(imdsEval.Labels); % number of images in dataset
numClasses = net.OutputSize; % # of classes in dataset

% Adversarial attack (L_inf attack)
% One way to define it is using original image +- disturbance (L_inf epsilon)
ones_ = ones([28 28], 'single'); % size of image
epsilon = 1; % pixel values (images are not normalized, they get normalized in the ImageInputLayer)

%% Main computation

% to save results (robustness and time)
results = zeros(N,2);
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

% Iterate trhough all images
for i=1:N

    % Load image in dataset
    [img, fileInfo] = readimage(imdsEval,i);
    target = single(fileInfo.Label); % label = 0 (index 1 for our network)
    img = single(img); % change precision

    % Adversarial attack
    disturbance = epsilon * ones_;
    % Ensure the values are within the valid range for pixels ([0 255])
    lb_min = zeros(size(img)); % minimum allowed values for lower bound is 0
    ub_max = 255*ones(size(img)); % maximum allowed values for upper bound is 255
    lb_clip = max((img-disturbance),lb_min);
    ub_clip = min((img+disturbance), ub_max);
    IS = ImageStar(lb_clip, ub_clip); % this is the input set we will use
    
    % Let's evaluate the image and the lower and upper bounds to ensure these
    % are correctly classified

    if ~mod(i,50)
        disp("Verifying image "+string(i)+" out of "+string(N)+" in the dataset...");
    end

    % Begin tracking time after input set is created
    t = tic;

    % Evaluate input image
    Y_outputs = net.evaluate(img); 
    [~, yPred] = max(Y_outputs);
    
    % Evaluate lower and upper bounds
    LB_outputs = net.evaluate(lb_clip);
    [~, LB_Pred] = max(LB_outputs); 
    UB_outputs = net.evaluate(ub_clip);
    [~, UB_Pred] = max(UB_outputs);

    % Check if outputs are violating robustness property
    if any([yPred, LB_Pred, UB_Pred] ~= target)
        results(i,1) = 0;
        results(i,2) = toc(t);
        % if counterexample found, no need to do any reachability analysis
        continue;
    end
    
    % If not, we verify the robustness using reachability analysis
    %  - Use the NN.verify_robustness function

    % A common approach would be to use some refinement approach like
    %   - Try first with faster approx method, if not robust, compute the
    %   exact reachability analysis

    % For the purpose of this tutorial, we are only going to do the
    % approximate method

    % Verification
    results(i,1) = net.verify_robustness(IS, reachOptions, target);
    results(i,2) = toc(t);

end

% Get summary results
N = length(results);
rob = sum(results(:,1) == 1);
not_rob = sum(results(:,1) == 0);
unk = sum(results(:,1) == 2);
totalTime = sum(results(:,2));
avgTime = totalTime/N;

% Print results to screen
disp("======= ROBUSTNESS RESULTS ==========")
disp(" ");
disp("Number of robust images = "+string(rob)+ ", equivalent to " + string(100*rob/N) + "% of the dataset.");
disp("Number of not robust images = " +string(not_rob)+ ", equivalent to " + string(100*not_rob/N) + "% of the dataset.")
disp("Number of unknown images = "+string(unk)+ ", equivalent to " + string(100*unk/N) + "% of the dataset.");
disp(" ");
disp("It took a total of "+string(totalTime) + " seconds to compute the verification results, an average of "+string(avgTime)+" seconds per image");

% Save results
save('results_verify_allTest.mat', 'results');