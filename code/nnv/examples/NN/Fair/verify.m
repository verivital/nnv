%% Fairness Verification of Adult Classification Model (NN)

% % Load network 
% ac_model = load('adult_model_fc.mat');
% 
% % Create NNV model 
% net = matlab2nnv(ac_model.net);

% 50_100_50 model
onnx_model_1 = fullfile('fair_model_ffnn_pytorch_50_100_50.onnx');

% Load the ONNX file as DAGNetwork
netONNX = importONNXNetwork(onnx_model_1, 'OutputLayerType', 'classification', 'InputDataFormats', {'BC'});

% Convert the DAGNetwork to NNV format
net = matlab2nnv(netONNX);
 
% Jimmy Rigged Fix: manually edit ouput size
net.OutputSize = 2;

% Load data
load("adult_fair_data.mat", 'X', 'y');
X_test_loaded = permute(X, [2, 1]); % change this for matlab expected format
y_test_loaded = y+1;  % update labels

% Normalize features in X_test_loaded
min_values = min(X_test_loaded, [], 2);
max_values = max(X_test_loaded, [], 2);

% Ensure no division by zero for constant features
variableFeatures = max_values - min_values > 0;
min_values(~variableFeatures) = 0; % Avoids changing constant features
max_values(~variableFeatures) = 1; % Avoids division by zero 

% Apply normalization
X_test_loaded(:, variableFeatures) = (X_test_loaded(:, variableFeatures) - min_values(variableFeatures)) ./ (max_values(variableFeatures) - min_values(variableFeatures));

% Count total observations
total_obs = size(X_test_loaded, 2);
% disp(['There are total ', num2str(total_obs), ' observations']);

% Number of observations we want to test
numObs = 200;

%% Verification

% to save results (robustness and time)
results = zeros(numObs,2);

% First, we define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approxiate method

nR = 67; % ---> just chosen arbitrarily

% ADJUST epsilon value here
% epsilon = [0.0001,0.01,0.25,0.5,1.0];
epsilon = [0.01,0.75];

%
% Set up results
%
nE = 1; %% will need to update later
res = zeros(numObs,nE); % robust result, changed from total_obs to 123
time = zeros(numObs,nE); % computation time, changed from total_obs to 123
met = repmat("relax", [numObs, nE]); % method used to compute result, changed from total_obs to 123

% Randomly select 20 observations
% rng(500); % Set a seed for reproducibility
rand_indices = randsample(total_obs, numObs);

for e=1:length(epsilon) 
% Iterate through all observations
for i=1:numObs
    idx = rand_indices(i);
    [IS, xRand] = L_inf_attack(X_test_loaded(:, idx), epsilon(e), nR, min_values, max_values);
     %
     % Try falsification, then relax star, if unknown, try approx-star
     %
     t = tic;
     for xR=1:length(nR+3)
        im = xRand(:, xR);
        predictedLabels = net.evaluate(im);
        [~, Pred] = max(predictedLabels);
        if Pred == y_test_loaded(i) % the current label of the orignal image
            res(i,e) = 0; % counterexample found
            time(i,e) = toc(t);
            met(i,e) = "counterexample";
            continue;
        end
     end

     try
        % relax star
        temp = net.verify_robustness(IS,reachOptions,y_test_loaded(i));
        % approx reachability if relax star unknown
        if temp ~= 1 && temp ~= 0
            reachOptions = struct;
            reachOptions.reachMethod = 'exact-star';
            temp = net.verify_robustness(IS,reachOptions,y_test_loaded(i));
            met(i,e) = 'exact';
        end
      catch ME
        met(i,e) = ME.message;
        temp = -1;
     end
      res(i,e) = temp; % robust result
      time(i, e) = toc(t); % store computation time

           
      % reset reachOptions
      reachOptions.reachMethod = 'relax-star-area';
      reachOptions.relaxFactor = 0.5;
end
% Get summary results
N = numObs;
rob = sum(res(:,e)==1);
not_rob = sum(res(:,e) == 0);
unk = sum(res(:,e) == 2);
totalTime = sum(time(:,e));
avgTime = totalTime/N;

% Print results to screen
disp("======= ROBUSTNESS RESULTS e: "+string(epsilon(e))+" ==========")
disp(" ");
disp("Number of robust samples = "+string(rob)+ ", equivalent to " + string(100*rob/N) + "% of the samples.");
disp("Number of not robust samples = " +string(not_rob)+ ", equivalent to " + string(100*not_rob/N) + "% of the samples.")
disp("Number of unknown samples = "+string(unk)+ ", equivalent to " + string(100*unk/N) + "% of the samples.");
disp(" ");
disp("It took a total of "+string(totalTime) + " seconds to compute the verification results, an average of "+string(avgTime)+" seconds per sample");

save("robustness_results_ex1","res","time","epsilon","met");
end

%% Helper Function
% Adjusted for fairness check -> only apply perturbation to desired
% sensitive feature.
function [IS, xRand] = L_inf_attack(x, epsilon, nR, min_values, max_values)
    % Apply a L-infinity attack of value epsilon on selected features of input image x
    % where VariableFeatures is true.
    % Return an ImageStar (IS) and random images from initial set
    SampleSize = size(x);

    % Initialize disturbance to only affect specified sensitive rows
    disturbance = zeros(SampleSize, "like", x);
    row_indices = [1, 9]; % Define which rows to perturb
    
    % Check if the row indices are within the valid range
    if max(row_indices) <= size(x, 1)  % Ensure the highest index doesn't exceed the number of rows
        disturbance(row_indices, :) = epsilon;  % Apply epsilon to specified rows across all columns
    else
        error('The input data does not have enough rows.');
    end


    % Calculate disturbed lower and upper bounds considering min and max values
    lb = max(x - disturbance, min_values);
    ub = min(x + disturbance, max_values);

    IS = ImageStar(single(lb), single(ub)); % default: single (assume onnx input models)

    % Create random images from initial set
    % Adjusted reshaping according to your specific needs
    lb = reshape(lb, [13,1]);  % Update the reshape parameters as per your actual data dimension
    ub = reshape(ub, [13,1]);
    xB = Box(single(lb), single(ub));
    xRand = xB.sample(nR);
    xRand = reshape(xRand,[13,nR]);
    xRand(:,nR+1) = x; % add original image
    xRand(:,nR+2) = IS.im_lb; % add lower bound image
    xRand(:,nR+3) = IS.im_ub; % add upper bound image
end

