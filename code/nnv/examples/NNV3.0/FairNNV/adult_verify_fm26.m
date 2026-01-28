%% Exact Fairness Verification of Adult Classification Model (NN)
% FM26 Tool Paper Results Script
% Generates results for: (1) Counterfactual fairness table
%                        (2) Individual fairness stacked bar charts
%                        (3) Comprehensive timing table
%
% This script can be run standalone or called from run_fm26_fairnnv.m
% When run standalone, default paths are used.
% When called from runner, paths are set via config struct.

% Suppress warnings
warning('off', 'nnet_cnn_onnx:onnx:WarnAPIDeprecation');
warning('off', 'nnet_cnn_onnx:onnx:FillingInClassNames');

%% Load data into NNV
warning('on', 'verbose')

%% Setup
% Check if config exists (set by runner script), otherwise use defaults
if ~exist('config', 'var')
    % Default configuration for standalone execution
    config.modelsDir = './adult_onnx';
    config.dataDir = './data';
    config.outputDir = './fm26_fairnnv_results';
    config.dataFile = 'adult_data.mat';
    config.modelList = {'AC-1', 'AC-3'};
    config.numObs = 100;
    config.randomSeed = 500;
    config.timeout = 600;
    config.epsilon_counterfactual = [0.0];
    config.epsilon_individual = [0.01, 0.02, 0.03, 0.05, 0.07, 0.1];

    % Clear workspace except config when running standalone
    clearvars -except config;
end

modelDir = config.modelsDir;
onnxFiles = dir(fullfile(modelDir, '*.onnx'));  % List all .onnx files

% Load data
dataFilePath = fullfile(config.dataDir, config.dataFile);
load(dataFilePath, 'X', 'y');

% Create results directory if it doesn't exist
resultsDir = config.outputDir;
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% Initialize results storage
results_counterfactual = {};  % For counterfactual fairness (epsilon = 0)
results_individual = {};      % For individual fairness (epsilon > 0)
results_timing = {};          % For comprehensive timing table

% List of models to process
modelList = config.modelList;

% Epsilon values
% 0.0 -> counterfactual fairness (flips sensitive attribute)
% >0.0 -> individual fairness (flips SA w/ perturbation of numerical features)
epsilon_counterfactual = config.epsilon_counterfactual;
epsilon_individual = config.epsilon_individual;
epsilon = [epsilon_counterfactual, epsilon_individual];  % Combined for processing

% Number of observations to test
numObs = config.numObs;

%% Loop through each model
for k = 1:length(onnxFiles)
    [~, modelName, ~] = fileparts(onnxFiles(k).name);
    if any(strcmp(modelName, modelList))

    clear net netONNX outputSet IS R

    onnx_model_path = fullfile(onnxFiles(k).folder, onnxFiles(k).name);

    % Load the ONNX file as DAGNetwork
    netONNX = importONNXNetwork(onnx_model_path, 'OutputLayerType', 'classification', 'InputDataFormats', {'BC'});

    % Convert the DAGNetwork to NNV format
    net = matlab2nnv(netONNX);

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

    % Test accuracy --> verify matches with python
    total_corr= 0;
    for i=1:total_obs
        im = X_test_loaded(:, i);
        predictedLabels = net.evaluate(im);
        [~, Pred] = min(predictedLabels);
        TrueLabel = y_test_loaded(i);
        if Pred == TrueLabel
            total_corr = total_corr + 1;
        end
    end
    disp("Model: " + modelName);
    disp("Accuracy of Model: "+string(total_corr/total_obs));

    %% Verification

    % First, we define the reachability options
    reachOptions = struct; % initialize
    reachOptions.reachMethod = 'exact-star';

    % Set up results
    nE = length(epsilon);
    res = zeros(numObs,nE); % robust result
    time = zeros(numObs,nE); % computation time
    met = repmat("exact", [numObs, nE]); % method used to compute result

    % Randomly select observations
    rng(config.randomSeed); % Set a seed for reproducibility
    rand_indices = randsample(total_obs, numObs);

    for e=1:length(epsilon)
        clear R outputSet IS temp t
        % Reset the timeout flag
        assignin('base', 'timeoutOccurred', false);

        % Create and configure the timer
        verificationTimer = timer;
        verificationTimer.StartDelay = config.timeout;  % Set timer (default: 10 minutes)
        verificationTimer.TimerFcn = @(myTimerObj, thisEvent) ...
        assignin('base', 'timeoutOccurred', true);
        start(verificationTimer);  % Start the timer


        for i=1:numObs
            idx = rand_indices(i);
            IS = perturbationIF(X_test_loaded(:, idx), epsilon(e), min_values, max_values);

            t = tic;  % Start timing the verification for each sample
            outputSet = net.reach(IS,reachOptions); % Generate output set
            target = y_test_loaded(idx);

            R = Star;
             % Process set
            if ~isa(outputSet, "Star")
                nr = length(outputSet);
                R(nr) = Star;
                for s=1:nr
                    R(s) = outputSet(s).toStar;
                end
            else
                R = outputSet;
            end

            % Process fairness specification
            target = net.robustness_set(target, 'min');

            % Verify fairness
            temp = 1;
            for s = 1:length(R)
                if verify_specification(R(s), target) ~= 1
                    temp = 0;
                    break
                end
            end

            met(i,e) = 'exact';

            res(i,e) = temp; % robust result

            time(i,e) = toc(t); % store computation time

            % Check for timeout flag
            if evalin('base', 'timeoutOccurred')
                disp(['Timeout reached for epsilon = ', num2str(epsilon(e)), ': stopping verification for this epsilon.']);
                res(i+1:end,e) = 2; % Mark remaining as unknown
                break; % Exit the inner loop after timeout
            end
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
        disp("======= FAIRNESS RESULTS e: "+string(epsilon(e))+" ==========")
        disp(" ");
        disp("Number of fair samples = "+string(rob)+ ", equivalent to " + string(100*rob/N) + "% of the samples.");
        disp("Number of non-fair samples = " +string(not_rob)+ ", equivalent to " + string(100*not_rob/N) + "% of the samples.")
        disp("Number of unknown samples = "+string(unk)+ ", equivalent to " + string(100*unk/N) + "% of the samples.");
        disp(" ");
        disp("It took a total of "+string(totalTime) + " seconds to compute the verification results, an average of "+string(avgTime)+" seconds per sample");

        % Collect results based on epsilon type
        if epsilon(e) == 0.0
            % Counterfactual fairness results
            results_counterfactual{end+1} = {modelName, 100 * rob / N, 100 * not_rob / N};
        else
            % Individual fairness results (for stacked bar chart)
            results_individual{end+1} = {modelName, epsilon(e), 100 * rob / N, 100 * not_rob / N, 100 * unk / N};
        end

        % Timing results (all epsilon values)
        results_timing{end+1} = {modelName, epsilon(e), totalTime, avgTime};
    end
    end
end

%% Save results to CSV files
% Get the current timestamp using datetime
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
timestampStr = char(timestamp);

% --- Save Counterfactual Fairness Results ---
% For Table: Counterfactual Fairness (epsilon = 0)
csv_counterfactual = fullfile(resultsDir, ['fm26_counterfactual_', timestampStr, '.csv']);
fid = fopen(csv_counterfactual, 'w');
fprintf(fid, 'Model,FairPercent,UnfairPercent\n');
for r = 1:length(results_counterfactual)
    fprintf(fid, '%s,%f,%f\n', results_counterfactual{r}{1}, results_counterfactual{r}{2}, results_counterfactual{r}{3});
end
fclose(fid);
disp(['Counterfactual results saved to ', csv_counterfactual]);

% --- Save Individual Fairness Results ---
% For Stacked Bar Charts (epsilon > 0)
csv_individual = fullfile(resultsDir, ['fm26_individual_', timestampStr, '.csv']);
fid = fopen(csv_individual, 'w');
fprintf(fid, 'Model,Epsilon,FairPercent,UnfairPercent,UnknownPercent\n');
for r = 1:length(results_individual)
    fprintf(fid, '%s,%f,%f,%f,%f\n', results_individual{r}{1}, results_individual{r}{2}, results_individual{r}{3}, results_individual{r}{4}, results_individual{r}{5});
end
fclose(fid);
disp(['Individual fairness results saved to ', csv_individual]);

% --- Save Comprehensive Timing Table ---
csv_timing = fullfile(resultsDir, ['fm26_timing_', timestampStr, '.csv']);
fid = fopen(csv_timing, 'w');
fprintf(fid, 'Model,Epsilon,TotalTime,AvgTimePerSample\n');
for r = 1:length(results_timing)
    fprintf(fid, '%s,%f,%f,%f\n', results_timing{r}{1}, results_timing{r}{2}, results_timing{r}{3}, results_timing{r}{4});
end
fclose(fid);
disp(['Timing results saved to ', csv_timing]);

disp(" ");
disp("======= FM26 VERIFICATION COMPLETE ==========");
disp("Generated files:");
disp("  1. " + csv_counterfactual + " (for counterfactual fairness table)");
disp("  2. " + csv_individual + " (for individual fairness stacked bar charts)");
disp("  3. " + csv_timing + " (for comprehensive timing table)");

%% Helper Function
% Apply perturbation (individual fairness) to sample
function IS = perturbationIF(x, epsilon, min_values, max_values)
    % Applies perturbations on selected features of input sample x
    % Return an ImageStar (IS) and random images from initial set
    SampleSize = size(x);

    disturbance = zeros(SampleSize, "like", x);
    sensitive_rows = [9];
    nonsensitive_rows = [1,10,11,12];

    % Flip the sensitive attribute
    if epsilon ~= -1
        if x(sensitive_rows) == 1
            x(sensitive_rows) = 0;
        else
            x(sensitive_rows) = 1;
        end
    end

    % Apply epsilon perturbation to non-sensitive numerical features
    for i = 1:length(nonsensitive_rows)
        if epsilon ~= -1
            if nonsensitive_rows(i) <= size(x, 1)
                disturbance(nonsensitive_rows(i), :) = epsilon;
            else
                error('The input data does not have enough rows.');
            end
        end
    end

    % Calculate disturbed lower and upper bounds considering min and max values
    lb = max(x - disturbance, min_values);
    ub = min(x + disturbance, max_values);
    IS = Star(single(lb), single(ub)); % default: single (assume onnx input models)

end
