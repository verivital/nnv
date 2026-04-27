%% FairNNV - Main Runner Script
% Runs the complete FairNNV verification pipeline on the Adult-Income
% classifier (counterfactual + individual fairness).
%
% USAGE:
%   matlab -batch "run_fairnnv"
%
% OUTPUTS (under FairNNV/results/<timestamp>/):
%   - CSV files with verification results
%   - PNG/PDF figures
%   - LaTeX tables (counterfactual, timing)
%
% REQUIREMENTS:
%   - NNV toolbox available at <repo>/code/nnv (resolved automatically below)
%   - ONNX models in FairNNV/models/ (AC-1.onnx, AC-3.onnx)
%   - adult_data.mat in FairNNV/data/

%% ================== CONFIGURATION ==================
% Defaults below are applied only for fields the caller has not already
% set in a pre-populated `config` struct (e.g., from run_all.sh --smoke).

if ~exist('config', 'var'); config = struct(); end
scriptDir = fileparts(mfilename('fullpath'));

% Path to NNV root folder (three levels up from FairNNV/).
if ~isfield(config, 'nnvDir');     config.nnvDir    = fullfile(scriptDir, '..', '..', '..'); end

% Bundled ONNX models and data live next to this script.
if ~isfield(config, 'modelsDir');  config.modelsDir = fullfile(scriptDir, 'models'); end
if ~isfield(config, 'dataDir');    config.dataDir   = fullfile(scriptDir, 'data'); end

% Timestamped output directory: results/<timestamp>/
if ~isfield(config, 'outputDir')
    ts = char(datetime('now', 'Format', 'yyMMdd-HHmmss'));
    config.outputDir = fullfile(scriptDir, 'results', ts);
end

% Data file name
if ~isfield(config, 'dataFile');   config.dataFile  = 'adult_data.mat'; end

% Models to verify (must match filenames without .onnx extension)
% AC-1: 13→16→8→2 (2 hidden, narrow)
% AC-3: 13→50→2 (1 hidden, shallow)
% AC-4: 13→100→100→2 (2 hidden, wide)
if ~isfield(config, 'modelList');  config.modelList = {'AC-1', 'AC-3'}; end

% Number of samples to test (default: 100)
if ~isfield(config, 'numObs');     config.numObs    = 100; end

% Random seed for reproducibility
if ~isfield(config, 'randomSeed'); config.randomSeed = 500; end

% Timeout per epsilon value in seconds (default: 600 = 10 minutes)
if ~isfield(config, 'timeout');    config.timeout   = 600; end

% Epsilon values for verification
% 0.0 = counterfactual fairness (flip sensitive attribute only)
% >0.0 = individual fairness (flip SA + perturb numerical features)
if ~isfield(config, 'epsilon_counterfactual'); config.epsilon_counterfactual = [0.0]; end
if ~isfield(config, 'epsilon_individual');     config.epsilon_individual = [0.01, 0.02, 0.03, 0.05, 0.07, 0.1]; end

% Figure export formats (set to false to skip)
if ~isfield(config, 'savePNG');    config.savePNG = true; end
if ~isfield(config, 'savePDF');    config.savePDF = true; end

%% ================== END CONFIGURATION ==================

%% Initialize NNV
disp("======= FairNNV Pipeline ==========");
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

run('adult_verify.m');

disp(" ");
disp("Verification complete.");
disp(" ");

%% Run Plotting
disp("======= STEP 2: Generating Figures ==========");
disp(" ");

run('plot_results.m');

disp(" ");
disp("======= FairNNV Pipeline Complete ==========");
disp(" ");
disp("All results saved to: " + config.outputDir);
