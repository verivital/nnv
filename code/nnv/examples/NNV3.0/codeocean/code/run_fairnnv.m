function run_fairnnv()
%% FairNNV Test Runner for CodeOcean
% Runs fairness verification on Adult Census dataset
%
% Data required in /data/:
%   - ICAIF24/adult_mat/AC-1.mat, AC-3.mat (pre-converted NNV networks)
%   - ICAIF24/data/adult_data.mat
%
% Output in /results/FairNNV/:
%   - fm26_counterfactual_*.csv
%   - fm26_individual_*.csv
%   - fm26_timing_*.csv

disp('Starting FairNNV verification...');

%% Suppress warnings related to loading pre-converted networks
% These warnings occur because the .mat files contain serialized MATLAB DL objects
% that reference classes not available without the ONNX support package.
% The warnings don't affect execution since we're using pre-converted NNV networks.
warning('off', 'all');  % Suppress all warnings during this function

%% Configuration for CodeOcean
config.modelsDir = '/data/ICAIF24/adult_mat';  % Use pre-converted .mat files
config.dataDir = '/data/ICAIF24/data';
config.outputDir = '/results/FairNNV';
config.dataFile = 'adult_data.mat';
config.modelList = {'AC-1', 'AC-3'};
config.numObs = 100;           % Number of samples to test
config.randomSeed = 500;       % For reproducibility
config.timeout = 600;          % 10 minutes per epsilon
config.epsilon_counterfactual = [0.0];
config.epsilon_individual = [0.01, 0.02, 0.03, 0.05, 0.07, 0.1];

%% Validate paths
if ~exist(config.modelsDir, 'dir')
    error('Models directory not found: %s\nPlease upload adult_mat folder to /data/ICAIF24/', config.modelsDir);
end

if ~exist(fullfile(config.dataDir, config.dataFile), 'file')
    error('Data file not found: %s\nPlease upload adult_data.mat to /data/ICAIF24/data/', fullfile(config.dataDir, config.dataFile));
end

% Create output directory
if ~exist(config.outputDir, 'dir')
    mkdir(config.outputDir);
end

%% Load data
disp('Loading data...');
modelDir = config.modelsDir;
matFiles = dir(fullfile(modelDir, '*.mat'));

dataFilePath = fullfile(config.dataDir, config.dataFile);
load(dataFilePath, 'X', 'y');

resultsDir = config.outputDir;
modelList = config.modelList;
epsilon_counterfactual = config.epsilon_counterfactual;
epsilon_individual = config.epsilon_individual;
epsilon = [epsilon_counterfactual, epsilon_individual];
numObs = config.numObs;

% Initialize results storage
results_counterfactual = {};
results_individual = {};
results_timing = {};

%% Loop through each model
for k = 1:length(matFiles)
    [~, modelName, ~] = fileparts(matFiles(k).name);
    if any(strcmp(modelName, modelList))

        disp(['Processing model: ' modelName]);

        matFilePath = fullfile(matFiles(k).folder, matFiles(k).name);
        disp(['  Loading from: ' matFilePath]);

        % Load pre-converted NNV network
        try
            loaded = load(matFilePath, 'net');
            net = loaded.net;
            disp(['  Loaded pre-converted network successfully']);
        catch loadErr
            disp(['  Error loading network: ' loadErr.message]);
            disp('  Attempting alternative load method...');
            % Try loading without specifying variable name
            loaded = load(matFilePath);
            if isfield(loaded, 'net')
                net = loaded.net;
                disp('  Alternative load succeeded');
            else
                rethrow(loadErr);
            end
        end

        % Prepare data
        X_test_loaded = permute(X, [2, 1]);
        y_test_loaded = y + 1;

        % Normalize features
        min_values = min(X_test_loaded, [], 2);
        max_values = max(X_test_loaded, [], 2);
        variableFeatures = max_values - min_values > 0;
        min_values(~variableFeatures) = 0;
        max_values(~variableFeatures) = 1;
        X_test_loaded = (X_test_loaded - min_values) ./ (max_values - min_values);

        total_obs = size(X_test_loaded, 2);

        % Test accuracy
        total_corr = 0;
        for i = 1:total_obs
            im = X_test_loaded(:, i);
            predictedLabels = net.evaluate(im);
            [~, Pred] = min(predictedLabels);
            TrueLabel = y_test_loaded(i);
            if Pred == TrueLabel
                total_corr = total_corr + 1;
            end
        end
        disp(['Accuracy of ' modelName ': ' num2str(total_corr/total_obs)]);

        %% Verification
        reachOptions = struct;
        reachOptions.reachMethod = 'exact-star';

        nE = length(epsilon);
        res = zeros(numObs, nE);
        time = zeros(numObs, nE);

        rng(config.randomSeed);
        rand_indices = randsample(total_obs, numObs);

        for e = 1:length(epsilon)
            assignin('base', 'timeoutOccurred', false);

            verificationTimer = timer;
            verificationTimer.StartDelay = config.timeout;
            verificationTimer.TimerFcn = @(myTimerObj, thisEvent) assignin('base', 'timeoutOccurred', true);
            start(verificationTimer);

            for i = 1:numObs
                idx = rand_indices(i);
                IS = perturbationIF(X_test_loaded(:, idx), epsilon(e), min_values, max_values);

                t = tic;
                outputSet = net.reach(IS, reachOptions);
                target = y_test_loaded(idx);

                R = Star;
                if ~isa(outputSet, "Star")
                    nr = length(outputSet);
                    R(nr) = Star;
                    for s = 1:nr
                        R(s) = outputSet(s).toStar;
                    end
                else
                    R = outputSet;
                end

                target = net.robustness_set(target, 'min');

                temp = 1;
                for s = 1:length(R)
                    if verify_specification(R(s), target) ~= 1
                        temp = 0;
                        break
                    end
                end

                res(i, e) = temp;
                time(i, e) = toc(t);

                if evalin('base', 'timeoutOccurred')
                    disp(['Timeout reached for epsilon = ' num2str(epsilon(e))]);
                    res(i+1:end, e) = 2;
                    break;
                end
            end

            stop(verificationTimer);
            delete(verificationTimer);

            % Summary
            N = numObs;
            rob = sum(res(:, e) == 1);
            not_rob = sum(res(:, e) == 0);
            unk = sum(res(:, e) == 2);
            totalTime = sum(time(:, e));
            avgTime = totalTime / N;

            disp(['Epsilon=' num2str(epsilon(e)) ': Fair=' num2str(rob) ' Unfair=' num2str(not_rob) ' Unknown=' num2str(unk)]);

            if epsilon(e) == 0.0
                results_counterfactual{end+1} = {modelName, 100 * rob / N, 100 * not_rob / N};
            else
                results_individual{end+1} = {modelName, epsilon(e), 100 * rob / N, 100 * not_rob / N, 100 * unk / N};
            end

            results_timing{end+1} = {modelName, epsilon(e), totalTime, avgTime};
        end
    end
end

%% Save results
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
timestampStr = char(timestamp);

% Counterfactual
csv_counterfactual = fullfile(resultsDir, ['fm26_counterfactual_' timestampStr '.csv']);
fid = fopen(csv_counterfactual, 'w');
fprintf(fid, 'Model,FairPercent,UnfairPercent\n');
for r = 1:length(results_counterfactual)
    fprintf(fid, '%s,%f,%f\n', results_counterfactual{r}{1}, results_counterfactual{r}{2}, results_counterfactual{r}{3});
end
fclose(fid);
disp(['Saved: ' csv_counterfactual]);

% Individual
csv_individual = fullfile(resultsDir, ['fm26_individual_' timestampStr '.csv']);
fid = fopen(csv_individual, 'w');
fprintf(fid, 'Model,Epsilon,FairPercent,UnfairPercent,UnknownPercent\n');
for r = 1:length(results_individual)
    fprintf(fid, '%s,%f,%f,%f,%f\n', results_individual{r}{1}, results_individual{r}{2}, results_individual{r}{3}, results_individual{r}{4}, results_individual{r}{5});
end
fclose(fid);
disp(['Saved: ' csv_individual]);

% Timing
csv_timing = fullfile(resultsDir, ['fm26_timing_' timestampStr '.csv']);
fid = fopen(csv_timing, 'w');
fprintf(fid, 'Model,Epsilon,TotalTime,AvgTimePerSample\n');
for r = 1:length(results_timing)
    fprintf(fid, '%s,%f,%f,%f\n', results_timing{r}{1}, results_timing{r}{2}, results_timing{r}{3}, results_timing{r}{4});
end
fclose(fid);
disp(['Saved: ' csv_timing]);

disp('FairNNV verification complete.');

end

%% Helper Function
function IS = perturbationIF(x, epsilon, min_values, max_values)
    SampleSize = size(x);
    disturbance = zeros(SampleSize, "like", x);
    sensitive_rows = [9];
    nonsensitive_rows = [1, 10, 11, 12];

    if epsilon ~= -1
        if x(sensitive_rows) == 1
            x(sensitive_rows) = 0;
        else
            x(sensitive_rows) = 1;
        end
    end

    for i = 1:length(nonsensitive_rows)
        if epsilon ~= -1
            if nonsensitive_rows(i) <= size(x, 1)
                disturbance(nonsensitive_rows(i), :) = epsilon;
            end
        end
    end

    lb = max(x - disturbance, min_values);
    ub = min(x + disturbance, max_values);
    IS = Star(single(lb), single(ub));
end
