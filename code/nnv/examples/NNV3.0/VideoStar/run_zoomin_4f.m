%% VideoStar ZoomIn-4f - Main Runner Script
% Runs a subset of the ZoomIn-4f video verification experiments.
%
% USAGE:
%   matlab -batch "cd code/nnv/examples/NNV3.0/VideoStar; run_zoomin_4f"
%
% OUTPUTS (under VideoStar/results/<timestamp>/):
%   - eps=<i>_255.csv: per-epsilon CSV with sample, result, time, method
%
% REQUIREMENTS:
%   - NNV toolbox available at <repo>/code/nnv (resolved automatically below)
%   - Bundled ONNX model in VideoStar/models/zoomin_4f.onnx
%   - Bundled data in VideoStar/data/ZoomIn/*.npy
%   - Bundled npy-matlab and vvn src in this folder

%% ================== CONFIGURATION ==================
% Defaults below are applied only for fields the caller has not already
% set in a pre-populated `config` struct (e.g., from run_all.sh --smoke).

if ~exist('config', 'var'); config = struct(); end
scriptDir = fileparts(mfilename('fullpath'));

% Path to NNV root folder (three levels up from VideoStar/).
if ~isfield(config, 'nnvDir');        config.nnvDir       = fullfile(scriptDir, '..', '..', '..'); end

% Bundled assets next to this script.
if ~isfield(config, 'modelsDir');     config.modelsDir    = fullfile(scriptDir, 'models'); end
if ~isfield(config, 'dataDir');       config.dataDir      = fullfile(scriptDir, 'data'); end
if ~isfield(config, 'npyMatlabDir');  config.npyMatlabDir = fullfile(scriptDir, 'npy-matlab'); end
if ~isfield(config, 'vvnSrcDir');     config.vvnSrcDir    = fullfile(scriptDir, 'src', 'vvn'); end

% Timestamped output: results/<yymmdd-HHmmss>/
if ~isfield(config, 'outputDir')
    ts = char(datetime('now', 'Format', 'yyMMdd-HHmmss'));
    config.outputDir = fullfile(scriptDir, 'results', ts);
end

% Verification settings
if ~isfield(config, 'dsType');        config.dsType       = 'zoom_in'; end       % 'zoom_in' or 'zoom_out'
if ~isfield(config, 'sampleLen');     config.sampleLen    = 4; end               % 4, 8, or 16
if ~isfield(config, 'verAlgorithm');  config.verAlgorithm = 'relax'; end         % 'relax' or 'approx'
if ~isfield(config, 'numClasses');    config.numClasses   = 10; end
if ~isfield(config, 'epsilon');       config.epsilon      = [1/255; 2/255; 3/255]; end
if ~isfield(config, 'timeout');       config.timeout      = 1800; end
if ~isfield(config, 'sampleIndices'); config.sampleIndices = 1:10; end

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

% GPU forward-compat (Blackwell / RTX 5090 under MATLAB R2024b).
% No-op on hosts without an NVIDIA GPU.
try
    parallel.gpu.enableCUDAForwardCompatibility(true);
catch
end

% Add npy-matlab to path for reading .npy files
if ~exist(config.npyMatlabDir, 'dir')
    error("npy-matlab directory not found: %s", config.npyMatlabDir);
end
addpath(config.npyMatlabDir);
disp("npy-matlab paths added successfully.");

% Add bundled vvn src directory to path (for verifyvideo.m)
addpath(config.vvnSrcDir);
disp("vvn verification functions added to path.");

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

%% Change to VideoStar directory so verifyvideo.m's relative paths
% (`data/ZoomIn/...`, `models/...`) resolve against this folder.
currentDir = pwd;
cd(scriptDir);
disp("Changed working directory to VideoStar for data loading.");
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
