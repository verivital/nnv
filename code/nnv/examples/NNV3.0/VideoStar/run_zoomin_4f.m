%% VideoStar ZoomIn-4f - Main Runner Script
% This script runs a subset of the ZoomIn-4f verification experiments.
%
% INSTRUCTIONS:
%   1. Set the paths in the CONFIGURATION section below
%   2. Run this script from the VideoStar directory
%   3. Results will be saved to the specified output folder
%
% OUTPUTS:
%   - CSV files with verification results for each epsilon value
%
% REQUIREMENTS:
%   - NNV toolbox must be installed
%   - ONNX models from FORMALISE2025/models/
%   - Data files from FORMALISE2025/data/
%   - npy-matlab for reading .npy files

%% ================== CONFIGURATION ==================
% Set these paths before running

% Detect local path automatically and resolve to absolute path
scriptDir = fileparts(mfilename('fullpath'));
oldDir = cd(fullfile(scriptDir, '..', '..', '..')); % Go up to code/nnv
config.nnvDir = pwd; % code/nnv (absolute)
cd(oldDir);

% Path to FORMALISE2025 submission directory
config.formalise2025Dir = fullfile(config.nnvDir, 'examples', 'Submission', 'FORMALISE2025');

% Path to folder containing ONNX models (e.g., zoomin_4f.onnx)
config.modelsDir = fullfile(config.formalise2025Dir, 'models');

% Path to folder containing data files (use local VideoStar data)
config.dataDir = fullfile(scriptDir, 'data');

% Path to npy-matlab for reading .npy files
config.npyMatlabDir = fullfile(config.formalise2025Dir, 'npy-matlab', 'npy-matlab');

% Path to output folder for results (will be created if it doesn't exist)
config.outputDir = '/tmp/results/VideoStar/ZoomIn/4';

% Verification settings
config.dsType = 'zoom_in';          % Dataset type: 'zoom_in' or 'zoom_out'
config.sampleLen = 4;               % Number of frames: 4, 8, or 16
config.verAlgorithm = 'relax';      % Verification algorithm: 'relax' or 'approx'

% Number of classes
config.numClasses = 10;

% Epsilon values for verification (perturbation sizes)
config.epsilon = [1/255; 2/255; 3/255];

% Timeout per sample in seconds (default: 1800 = 30 minutes)
config.timeout = 1800;

% Sample indices to verify (subset for quick testing)
% Using first 10 samples for a subset run
config.sampleIndices = 1:10;

%% ================== END CONFIGURATION ==================

%% Initialize NNV
disp("======= VideoStar ZoomIn-4f Verification ==========");
disp(" ");
disp("Initializing NNV toolbox...");

% Check if NNV directory exists
if ~exist(config.nnvDir, 'dir')
    error("NNV directory not found: %s", config.nnvDir);
end

% Add NNV to path
addpath(genpath(config.nnvDir));
disp("NNV paths added successfully.");

% Add npy-matlab to path for reading .npy files
if ~exist(config.npyMatlabDir, 'dir')
    error("npy-matlab directory not found: %s", config.npyMatlabDir);
end
addpath(config.npyMatlabDir);
disp("npy-matlab paths added successfully.");

% Add FORMALISE2025 src directory to path (for verifyvideo.m)
addpath(fullfile(config.formalise2025Dir, 'src', 'vvn'));
disp("FORMALISE2025 verification functions added to path.");

disp(" ");

%% Validate Configuration
disp("Validating configuration...");

% Check if models directory exists
if ~exist(config.modelsDir, 'dir')
    error("Models directory not found: %s", config.modelsDir);
end

% Check if data directory exists
if ~exist(config.dataDir, 'dir')
    error("Data directory not found: %s", config.dataDir);
end

% Check if model exists
modelName = sprintf("zoomin_%df.onnx", config.sampleLen);
modelPath = fullfile(config.modelsDir, modelName);
if ~exist(modelPath, 'file')
    error("Model not found: %s", modelPath);
end
disp("  Found model: " + modelName);

% Create output directory if it doesn't exist
if ~exist(config.outputDir, 'dir')
    mkdir(config.outputDir);
    disp("  Created output directory: " + config.outputDir);
else
    disp("  Output directory: " + config.outputDir);
end

disp(" ");
disp("Configuration validated successfully.");
disp(" ");

%% Change to FORMALISE2025 directory (required for data loading paths)
currentDir = pwd;
cd(config.formalise2025Dir);
disp("Changed to FORMALISE2025 directory for data loading.");
disp(" ");

%% Run Verification
disp("======= Running ZoomIn-4f Verification ==========");
disp(" ");
disp(sprintf("Dataset: %s", config.dsType));
disp(sprintf("Sample Length: %d frames", config.sampleLen));
disp(sprintf("Verification Algorithm: %s", config.verAlgorithm));
disp(sprintf("Number of samples: %d", length(config.sampleIndices)));
disp(sprintf("Epsilon values: [1/255, 2/255, 3/255]"));
disp(sprintf("Timeout: %d seconds", config.timeout));
disp(" ");

% Initialize results storage
numSamples = length(config.sampleIndices);
numEpsilons = length(config.epsilon);
results = cell(numSamples, numEpsilons);
times = cell(numSamples, numEpsilons);
methods = cell(numSamples, numEpsilons);

% Run verification for each sample and epsilon
for sampleIdx = 1:numSamples
    sampleNum = config.sampleIndices(sampleIdx);
    disp(sprintf("Processing sample %d (%d/%d)...", sampleNum, sampleIdx, numSamples));

    for epsIdx = 1:numEpsilons
        eps = config.epsilon(epsIdx);
        disp(sprintf("  Epsilon = %d/255...", epsIdx));

        try
            % Call the verification function from FORMALISE2025
            [res, t, met] = verifyvideo(config.dsType, config.sampleLen, ...
                                        config.verAlgorithm, sampleNum, epsIdx);

            results{sampleIdx, epsIdx} = res;
            times{sampleIdx, epsIdx} = t;
            methods{sampleIdx, epsIdx} = met;

            disp(sprintf("    Result: %d, Time: %.2f s", res, t));

        catch ME
            disp(sprintf("    Error: %s", ME.message));
            results{sampleIdx, epsIdx} = -1;
            times{sampleIdx, epsIdx} = -1;
            methods{sampleIdx, epsIdx} = ME.message;
        end
    end
end

%% Save Results
disp(" ");
disp("======= Saving Results ==========");

% Save results for each epsilon value
for epsIdx = 1:numEpsilons
    outputFile = fullfile(config.outputDir, sprintf('eps=%d_255.csv', epsIdx));

    % Open file and write header
    fid = fopen(outputFile, 'w');
    fprintf(fid, 'Sample Number,Result,Time,Method\n');

    % Write results
    for sampleIdx = 1:numSamples
        sampleNum = config.sampleIndices(sampleIdx);
        res = results{sampleIdx, epsIdx};
        t = times{sampleIdx, epsIdx};
        met = methods{sampleIdx, epsIdx};

        if isnumeric(t)
            fprintf(fid, '%d,%d,%.6f,%s\n', sampleNum, res, t, met);
        else
            fprintf(fid, '%d,%d,%s,%s\n', sampleNum, res, t, met);
        end
    end

    fclose(fid);
    disp(sprintf("  Saved: %s", outputFile));
end

%% Display Summary
disp(" ");
disp("======= Verification Summary ==========");
disp(" ");

for epsIdx = 1:numEpsilons
    verified = 0;
    unknown = 0;
    violated = 0;
    timeout = 0;
    totalTime = 0;
    validTimes = 0;

    for sampleIdx = 1:numSamples
        res = results{sampleIdx, epsIdx};
        t = times{sampleIdx, epsIdx};

        if res == 1
            verified = verified + 1;
        elseif res == 0
            violated = violated + 1;
        elseif res == 2
            unknown = unknown + 1;
        elseif res == 3
            timeout = timeout + 1;
        end

        if isnumeric(t) && t > 0
            totalTime = totalTime + t;
            validTimes = validTimes + 1;
        end
    end

    avgTime = totalTime / max(validTimes, 1);

    disp(sprintf("Epsilon = %d/255:", epsIdx));
    disp(sprintf("  Verified (robust): %d", verified));
    disp(sprintf("  Violated: %d", violated));
    disp(sprintf("  Unknown: %d", unknown));
    disp(sprintf("  Timeout: %d", timeout));
    disp(sprintf("  Average time: %.2f s", avgTime));
    disp(" ");
end

%% Cleanup
cd(currentDir);

disp("======= VideoStar ZoomIn-4f Complete ==========");
disp(" ");
disp(sprintf("Results saved to: %s", config.outputDir));
