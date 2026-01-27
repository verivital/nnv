%% FM26 FairNNV - Main Runner Script
% This script runs the complete FairNNV verification pipeline for FM26.
%
% INSTRUCTIONS:
%   1. Set the paths in the CONFIGURATION section below
%   2. Run this script
%   3. Results will be saved to the specified output folder
%
% OUTPUTS:
%   - CSV files with verification results
%   - PNG/PDF figures for the paper
%   - Timing table
%
% REQUIREMENTS:
%   - NNV toolbox must be installed
%   - ONNX models in the specified models folder
%   - Data file (adult_data.mat) in the specified data folder

%% ================== CONFIGURATION ==================
% Set these paths before running

% Path to NNV root folder (contains startup_nnv.m)
% This is used to add NNV to the MATLAB path
% Detect local path automatically and resolve to absolute path
scriptDir = fileparts(mfilename('fullpath'));
oldDir = cd(fullfile(scriptDir, '..', '..', '..')); % Go up to code/nnv/examples -> code/nnv
config.nnvDir = pwd; % code/nnv (absolute)
cd(oldDir);

% Path to folder containing ONNX models (e.g., AC-1.onnx, AC-3.onnx)
config.modelsDir = fullfile(config.nnvDir, 'examples', 'Submission', 'ICAIF24', 'adult_onnx');

% Path to folder containing data files (e.g., adult_data.mat)
config.dataDir = fullfile(config.nnvDir, 'examples', 'Submission', 'ICAIF24', 'data');

% Path to output folder for results (will be created if it doesn't exist)
config.outputDir = './fm26_fairnnv_results';

% Data file name
config.dataFile = 'adult_data.mat';

% Models to verify (must match filenames without .onnx extension)
% AC-1: 13→16→8→2 (2 hidden, narrow)
% AC-3: 13→50→2 (1 hidden, shallow)
% AC-4: 13→100→100→2 (2 hidden, wide)
config.modelList = {'AC-1', 'AC-3'};

% Number of samples to test (default: 100)
config.numObs = 100;

% Random seed for reproducibility
config.randomSeed = 500;

% Timeout per epsilon value in seconds (default: 600 = 10 minutes)
config.timeout = 600;

% Epsilon values for verification
% 0.0 = counterfactual fairness (flip sensitive attribute only)
% >0.0 = individual fairness (flip SA + perturb numerical features)
config.epsilon_counterfactual = [0.0];
config.epsilon_individual = [0.01, 0.02, 0.03, 0.05, 0.07, 0.1];

% Figure export formats (set to false to skip)
config.savePNG = true;
config.savePDF = true;

%% ================== END CONFIGURATION ==================

%% Initialize NNV
disp("======= FM26 FairNNV Pipeline ==========");
disp(" ");
disp("Initializing NNV toolbox...");

% Check if NNV directory exists
if ~exist(config.nnvDir, 'dir')
    error("NNV directory not found: %s", config.nnvDir);
end

% Add NNV to path (simplified - avoids tbxmanager permission issues)
addpath(genpath(config.nnvDir));
disp("NNV paths added successfully.");
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

% Check if data file exists
dataFilePath = fullfile(config.dataDir, config.dataFile);
if ~exist(dataFilePath, 'file')
    error("Data file not found: %s", dataFilePath);
end

% Check if at least one model exists
modelFound = false;
for i = 1:length(config.modelList)
    modelPath = fullfile(config.modelsDir, [config.modelList{i}, '.onnx']);
    if exist(modelPath, 'file')
        modelFound = true;
        disp("  Found model: " + config.modelList{i});
    else
        warning("Model not found: %s", modelPath);
    end
end
if ~modelFound
    error("No models found in: %s", config.modelsDir);
end

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

%% Run Verification
disp("======= STEP 1: Running Verification ==========");
disp(" ");

run('adult_verify_fm26.m');

disp(" ");
disp("Verification complete.");
disp(" ");

%% Run Plotting
disp("======= STEP 2: Generating Figures ==========");
disp(" ");

run('plot_fm26_results.m');

disp(" ");
disp("======= FM26 FairNNV Pipeline Complete ==========");
disp(" ");
disp("All results saved to: " + config.outputDir);
