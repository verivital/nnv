%% Fairness Verification of Adult Classification Model (NN)
% Comparison for the models used in Fairify

% Suppress warnings
warning('off', 'nnet_cnn_onnx:onnx:WarnAPIDeprecation');
warning('off', 'nnet_cnn_onnx:onnx:FillingInClassNames');

%% Load data into NNV
warning('on', 'verbose')

%% Setup
clear; clc;
modelDir = './adult_onnx';  % Directory containing ONNX models
onnxFiles = dir(fullfile(modelDir, '*.onnx'));  % List all .onnx files

load("adult_data.mat", 'X', 'y'); % load Adult data

%% Loop through each model
for k = 1:length(onnxFiles)

    onnx_model_path = fullfile(onnxFiles(k).folder, onnxFiles(k).name);

    % Load the ONNX file as DAGNetwork
    netONNX = importONNXNetwork(onnx_model_path, 'OutputLayerType', 'classification', 'InputDataFormats', {'BC'});

    % Convert the DAGNetwork to NNV format
    net = matlab2nnv(netONNX);
     
    % Jimmy Rigged Fix: manually edit ouput size
    net.OutputSize = 2;
    
    X_test_loaded = permute(X, [2, 1]);
    y_test_loaded = y+1;  % update labels
    
    % Normalize features in X_test_loaded
    min_values = min(X_test_loaded, [], 2);
    max_values = max(X_test_loaded, [], 2);
    
    % Ensure no division by zero for constant features
    variableFeatures = max_values - min_values > 0;
    min_values(~variableFeatures) = 0; % Avoids changing constant features
    max_values(~variableFeatures) = 1; % Avoids division by zero 

    % Normalizing X_test_loaded
    X_test_loaded = (X_test_loaded - min_values) ./ (max_values - min_values);

    % Count total observations
    total_obs = size(X_test_loaded, 2);
    % disp(['There are total ', num2str(total_obs), ' observations']);
   
    % Number of observations we want to test
    numObs = 100;
    
    %% Verification
    
    % To save results (robustness and time)
    results = zeros(numObs,2);
    
    % First, we define the reachability options
    reachOptions = struct; % initialize
    reachOptions.reachMethod = 'relax-star-area'; % using relax method
    reachOptions.relaxFactor = 0.5;
    
    nR = 50; % ---> just chosen arbitrarily
    
    % ADJUST epsilon values here
    epsilon = [0.0,0.001,0.01];

   
    % Set up results
    nE = 3; 
    res = zeros(numObs,nE); % robust result
    time = zeros(numObs,nE); % computation time
    met = repmat("relax", [numObs, nE]); % method used to compute result
 
    
    % Randomly select observations
    rng(500); % Set a seed for reproducibility
    rand_indices = randsample(total_obs, numObs);
    
    for e=1:length(epsilon)
        % Reset the timeout flag
        assignin('base', 'timeoutOccurred', false);

        % Create and configure the timer
        verificationTimer = timer;
        verificationTimer.StartDelay = 600;  % Set timer for 10 minutes
        verificationTimer.TimerFcn = @(myTimerObj, thisEvent) ...
        assignin('base', 'timeoutOccurred', true);
        start(verificationTimer);  % Start the timer
        
        % Iterate through observations
        for i=1:numObs
            idx = rand_indices(i);
            [IS, xRand] = perturbationIF(X_test_loaded(:, idx), epsilon(e), nR, min_values, max_values);
           
            skipTryCatch = false;  % Initialize the flag
            t = tic;  % Start timing the verification for each sample
            %
            % Try falsification, then relax star, if unknown, try exact-star
            %
            % Falsification
            for xR=1:length(nR+3)
                im = xRand(:, xR);
                predictedLabels = net.evaluate(im);
                [~, Pred] = min(predictedLabels);
                if Pred ~= y_test_loaded(idx)
                    res(i,e) = 0; % counterexample found
                    time(i,e) = toc(t);
                    met(i,e) = "counterexample";
                    skipTryCatch = true;  % Set the flag to skip try-catch block
                    continue;
                end
            end
            if ~skipTryCatch
                try
                    % Relax star
                    temp = net.verify_robustness(IS, reachOptions, y_test_loaded(idx));
                    % Exact star
                    if temp ~= 1 && temp ~= 0
                        reachOptions.reachMethod = 'exact-star';
                        temp = net.verify_robustness(IS, reachOptions, y_test_loaded(idx));
                        met(i,e) = 'exact';
                    end
                    res(i,e) = temp; % robust result
                catch ME
                    met(i,e) = ME.message;
                    temp = -1;
                end
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
        fprintf('Model: %s\n', onnxFiles(k).name);
        disp("======= ROBUSTNESS RESULTS e: "+string(epsilon(e))+" ==========")
        disp(" ");
        disp("Number of fair samples = "+string(rob)+ ", equivalent to " + string(100*rob/N) + "% of the samples.");
        disp("Number of non-fair samples = " +string(not_rob)+ ", equivalent to " + string(100*not_rob/N) + "% of the samples.")
        disp("Number of unknown samples = "+string(unk)+ ", equivalent to " + string(100*unk/N) + "% of the samples.");
        disp(" ");
        disp("It took a total of "+string(totalTime) + " seconds to compute the verification results, an average of "+string(avgTime)+" seconds per sample");
    end
end


%% Helper Function
% Apply perturbation (individual fairness) to sample
function [IS, xRand] = perturbationIF(x, epsilon, nR, min_values, max_values)
    % Applies perturbations on selected features of input sample x
    % Return an ImageStar (IS) and random images from initial set
    SampleSize = size(x);

    disturbance = zeros(SampleSize, "like", x);
    sensitive_rows = [9]; 
    nonsensitive_rows = [1,10,11,12];
    
    % Flip the sensitive attribute
    if x(sensitive_rows) == 1
        x(sensitive_rows) = 0;
    else
        x(sensitive_rows) = 1;
    end

    % Apply epsilon perturbation to non-sensitive numerical features
    for i = 1:length(nonsensitive_rows)
        if nonsensitive_rows(i) <= size(x, 1)
            disturbance(nonsensitive_rows(i), :) = epsilon;
        else
            error('The input data does not have enough rows.');
        end
    end

    % Calculate disturbed lower and upper bounds considering min and max values
    lb = max(x - disturbance, min_values);
    ub = min(x + disturbance, max_values);
    IS = ImageStar(single(lb), single(ub)); % default: single (assume onnx input models)

    % Create random samples from initial set
    % Adjusted reshaping according to specific needs
    lb = reshape(lb, [13,1]); 
    ub = reshape(ub, [13,1]);
    xB = Box(single(lb), single(ub));
    xRand = xB.sample(nR);
    xRand = reshape(xRand,[13,nR]);
    xRand(:,nR+1) = x; % add original image
    xRand(:,nR+2) = IS.im_lb; % add lower bound image
    xRand(:,nR+3) = IS.im_ub; % add upper bound image
end

