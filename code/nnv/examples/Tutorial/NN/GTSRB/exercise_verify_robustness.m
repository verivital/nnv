%% Robustness verification of a NN (L infinity adversarial attack)
%  if f(x) = y, then forall x' in X s.t. ||x - x'||_{\infty} <= eps,
%  then f(x') = y = f(x)

% Load network (gtsrb_model.mat)
model = ; % load the file containing the neural network (net -> SeriesNetwork)

% Create NNV model
net = ; % convert Series Network (from MATLAB) to NNV

% Load data
gtsrb_path = [nnvroot(), filesep, 'data', filesep, 'GTSRB', filesep];
imds = imageDatastore(gtsrb_path, 'IncludeSubfolders',true,'LabelSource','foldernames');

inputSize = [30 29];
imds.ReadFcn = @(loc)imresize(imread(loc),inputSize);

% Load one image in dataset
im_idx = ; % choose an index number
[img, fileInfo] = readimage(imds,im_idx);

% Visualize image 
figure;
imshow(img);

% Get image info
target = double(fileInfo.Label); % target label
img = double(img); % convert to double
    
% Create input set

% need to define lower and upper bounds based on an L_infinity perturbation
% of epsilon value = 1
% Set must be an Imagestar defined by upper and lower bounds
IS = ; % this is the input set we will use (ImageStar)


% Let's evaluate the image and the lower and upper bounds to ensure these
% are correctly classified

% Evaluate input image
t = tic;
Y_outputs = ; % complete function to evaluate image
[~, yPred] = max(Y_outputs); % (expected: y = target)

% Evaluate lower and upper bounds
LB_outputs = ; % complete function to evaluate lower bound
[~, LB_Pred] = max(LB_outputs); % (expected: y = target)
UB_outputs = ; % complete function to evaluate upper bounds
[~, UB_Pred] = max(UB_outputs); % (expected: y = target)
t_eval = toc(t);

if yPred == target && any([LB_Pred, UB_Pred] ~= target)
    disp("Neural Network is not robust!");
    disp("Counterexample found in "+string(t_eval)+" seconds");
end

% Now, we can do the verification process of this image w/ L_inf attack

% The easiest way to do it is using the verify_robustness function

% First, we need to define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approxiate method

% Robustness verification
t = tic;
res_approx = net.verify_robustness(); % complete the function 

if res_approx == 1
    disp("Neural network is verified to be robust!")
elseif any([yPred, LB_Pred, UB_Pred] ~= target)
    disp("Neural network is not robust!");
else
    disp("Unknown result")
end

toc(t);

