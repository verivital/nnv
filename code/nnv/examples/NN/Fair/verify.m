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

analyzeNetwork(netONNX)

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

% % Apply normalization
% X_test_loaded(:, variableFeatures) = (X_test_loaded(:, variableFeatures) - min_values(variableFeatures)) ./ (max_values(variableFeatures) - min_values(variableFeatures));

% Count total observations
total_obs = size(X_test_loaded, 2);
% disp(['There are total ', num2str(total_obs), ' observations']);

% Number of observations we want to test
numObs = 100;

%% Verification

% to save results (robustness and time)
results = zeros(numObs,2);

% First, we define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'relax-star-area'; % using approxiate method

nR = 75; % ---> just chosen arbitrarily

% ADJUST epsilon value here
epsilon = [0.0];
% epsilon = [0.001,0.01,0.25];

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
    % Reset the timeout flag
    assignin('base', 'timeoutOccurred', false);

    % Create and configure the timer
    verificationTimer = timer;
    verificationTimer.StartDelay = 500;  % Set timer for 200 seconds
    verificationTimer.TimerFcn = @(myTimerObj, thisEvent) ...
    assignin('base', 'timeoutOccurred', true);
    start(verificationTimer);  % Start the timer

    for i=1:numObs
        idx = rand_indices(i);
        [IS, xRand] = L_inf_attack(X_test_loaded(:, idx), epsilon(e), nR, min_values, max_values);
       
        t = tic;  % Start timing the verification for each sample

        try
            for xR=1:length(nR+3)
                im = xRand(:, xR);
                predictedLabels = net.evaluate(im);
                [~, Pred] = max(predictedLabels);
                if Pred == y_test_loaded(i)
                    res(i,e) = 0; % counterexample found
                    time(i,e) = toc(t);
                    met(i,e) = "counterexample";
                    continue;
                end
            end

            temp = net.verify_robustness(IS, reachOptions, y_test_loaded(i));
            if temp ~= 1 && temp ~= 0
                reachOptions.reachMethod = 'approx-star';
                temp = net.verify_robustness(IS, reachOptions, y_test_loaded(i));
                met(i,e) = 'approx';
            end
            res(i,e) = temp; % robust result
        catch ME
            met(i,e) = ME.message;
            temp = -1;
        end

        time(i,e) = toc(t); % store computation time

        % Check for timeout flag
        if evalin('base', 'timeoutOccurred')
            disp(['Timeout reached for epsilon = ', num2str(epsilon(e)), ': stopping verification for this epsilon.']);
            res(i+1:end,e) = 2; % Mark remaining as unknown
            break; % Exit the inner loop after timeout
        end

        % Reset reachOptions
        reachOptions.reachMethod = 'relax-star-area';
        reachOptions.relaxFactor = 0.5;
    end

    % Summary results, stopping, and deleting the timer should be outside the inner loop
    stop(verificationTimer);
    delete(verificationTimer);

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
    specific_rows = [9]; % Rows with specific perturbations
    specific_epsilons = [1]; % Corresponding specific perturbations

    % Apply specific perturbations to rows 8 and 9
    for i = 1:length(specific_rows)
        if specific_rows(i) <= size(x, 1)
            disturbance(specific_rows(i), :) = specific_epsilons(i);
        else
            error('The input data does not have enough rows.');
        end
    end
    
    % Apply general epsilon to all other rows except 8 and 9
    for row = 1:size(x, 1)
        if ~ismember(row, specific_rows)
            disturbance(row, :) = epsilon;
        end
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

