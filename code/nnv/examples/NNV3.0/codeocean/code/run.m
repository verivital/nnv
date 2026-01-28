%% NNV 3.0 Test Suite - CodeOcean Main Runner
% This script runs all five NNV 3.0 verification test categories:
% 1. FairNNV - Fairness verification
% 2. ProbVer - Probabilistic verification (CP-Star)
% 3. ModelStar - Weight perturbation verification
% 4. VideoStar - Video classification verification
% 5. GNNV - Graph Neural Network verification
%
% CodeOcean Paths:
%   /code/    - This script and test runners
%   /data/    - Input data files
%   /results/ - Output results

%% Setup
disp('========================================');
disp('   NNV 3.0 Verification Test Suite');
disp('   CodeOcean Capsule');
disp('========================================');
disp(' ');

% Add dependencies from /deps if exists (CodeOcean postInstall)
if exist('/deps', 'dir')
    addpath(genpath('/deps'));
    disp('Added /deps to path.');
end

%% Install ONNX support package at runtime (for ProbVer and VideoStar)
% Note: /code and /data are only available at runtime, not during Docker build
sppRoot = fullfile(matlabroot, 'SupportPackages');
sppInstalled = false;

% Check for sppFile.zip in /code/support_packages or /data/support_packages
sppZipPaths = {'/code/support_packages/sppFile.zip', '/data/support_packages/sppFile.zip'};
for i = 1:length(sppZipPaths)
    if exist(sppZipPaths{i}, 'file')
        disp(['Found support package: ' sppZipPaths{i}]);
        disp('Installing ONNX support package...');

        % Create SupportPackages directory
        if ~exist(sppRoot, 'dir')
            mkdir(sppRoot);
        end

        % Unzip the support package
        try
            unzip(sppZipPaths{i}, sppRoot);
            disp('Support package extracted successfully.');
            sppInstalled = true;
        catch ME
            disp(['Warning: Could not extract support package: ' ME.message]);
        end
        break;
    end
end

if ~sppInstalled
    disp('No sppFile.zip found - ProbVer and VideoStar require ONNX support package.');
end

% Set support package root if available
try
    if exist(sppRoot, 'dir')
        matlabshared.supportpkg.setSupportPackageRoot(sppRoot);
        addpath(genpath(sppRoot));
        disp(['Support package root: ' sppRoot]);
    end
catch ME
    disp(['Warning: Could not set support package root: ' ME.message]);
end

% Add /code to path (contains stub ONNX layer classes for FairNNV)
if exist('/code', 'dir')
    addpath('/code');
    disp('Added /code to path (includes ONNX stub classes).');
end

% Suppress message catalog warnings
warning('off', 'MATLAB:msgcat:CatalogNotFound');

% Add NNV to path - check multiple possible locations
nnvPath = '';
possiblePaths = {'/code/nnv', '/home/matlab/nnv/code/nnv', '/nnv/code/nnv', '/tmp/nnv/code/nnv'};
for i = 1:length(possiblePaths)
    if exist(possiblePaths{i}, 'dir')
        nnvPath = possiblePaths{i};
        break;
    end
end

% If NNV not found, clone it
if isempty(nnvPath)
    disp('NNV not found. Cloning repository...');
    nnvPath = '/tmp/nnv/code/nnv';
    [status, result] = system('git clone --depth 1 https://github.com/verivital/nnv.git /tmp/nnv');
    if status ~= 0
        error('Failed to clone NNV: %s', result);
    end
    disp('NNV cloned successfully.');
end

addpath(genpath(nnvPath));
disp(['NNV path added: ' nnvPath]);

% Add npy-matlab to path (for VideoStar) - check multiple locations
npyPath = '';
possibleNpyPaths = {'/code/npy-matlab/npy-matlab', '/home/matlab/npy-matlab/npy-matlab', '/npy-matlab/npy-matlab', '/tmp/npy-matlab/npy-matlab'};
for i = 1:length(possibleNpyPaths)
    if exist(possibleNpyPaths{i}, 'dir')
        npyPath = possibleNpyPaths{i};
        break;
    end
end

% If npy-matlab not found, clone it
if isempty(npyPath)
    disp('npy-matlab not found. Cloning repository...');
    npyPath = '/tmp/npy-matlab/npy-matlab';
    [status, result] = system('git clone --depth 1 https://github.com/kwikteam/npy-matlab.git /tmp/npy-matlab');
    if status ~= 0
        warning('Failed to clone npy-matlab: %s. VideoStar may not work.', result);
    else
        disp('npy-matlab cloned successfully.');
    end
end

if exist(npyPath, 'dir')
    addpath(npyPath);
    disp(['npy-matlab path added: ' npyPath]);
end

% Add YAML library to path (for ModelStar)
yamlPath = fullfile(nnvPath, 'engine', 'utils', 'yaml');
if exist(yamlPath, 'dir')
    addpath(yamlPath);
    disp('YAML library path added.');
end

% Create results directories
resultsDir = '/results';
mkdir(fullfile(resultsDir, 'FairNNV'));
mkdir(fullfile(resultsDir, 'ProbVer'));
mkdir(fullfile(resultsDir, 'ModelStar'));
mkdir(fullfile(resultsDir, 'VideoStar'));
mkdir(fullfile(resultsDir, 'GNNV'));
disp('Results directories created.');

%% Setup Python virtual environment for ProbVer (CP-Star verification)
% The cp_env.m function looks for .venv relative to nnvroot()
% nnvroot() returns the path 4 levels up from code/nnv/engine/utils/nnvroot.m
% For CodeOcean with NNV at /code/nnv, nnvroot() returns / (root)
% So we need to create the venv at the location where nnvroot()/.venv will find it
disp(' ');
disp('Setting up Python virtual environment for ProbVer...');

% Determine where cp_env.m will look for the venv
try
    nnv_root = nnvroot();
    expectedVenvPath = fullfile(nnv_root, '.venv');
    disp(['  NNV root: ' nnv_root]);
    disp(['  Expected venv path: ' expectedVenvPath]);
catch
    % If nnvroot fails, default to /code/.venv
    nnv_root = '/code';
    expectedVenvPath = '/code/.venv';
    disp(['  NNV root detection failed, using default: ' nnv_root]);
end

% Check multiple possible venv locations (in priority order)
venvLocations = {
    expectedVenvPath, ...              % Where cp_env.m expects it
    '/deps/probver_venv', ...          % Created by postInstall
    '/code/.venv', ...                 % Alternative location
    '/tmp/.venv'                       % Last resort
};

venvPath = '';
pythonExe = '';
for i = 1:length(venvLocations)
    testPath = venvLocations{i};
    testPython = fullfile(testPath, 'bin', 'python');
    if exist(testPython, 'file')
        venvPath = testPath;
        pythonExe = testPython;
        disp(['  Found existing venv: ' venvPath]);
        break;
    end
end

% If no venv found, create one
if isempty(venvPath)
    disp('  No existing venv found, creating new one...');

    % Try to create at expected location first
    for i = 1:length(venvLocations)
        targetPath = venvLocations{i};
        disp(['  Trying to create venv at: ' targetPath]);

        [status, result] = system(['python3 -m venv "' targetPath '" 2>&1']);
        if status == 0
            venvPath = targetPath;
            pythonExe = fullfile(venvPath, 'bin', 'python');
            disp(['  Created venv: ' venvPath]);
            break;
        else
            disp(['    Failed: ' strtrim(result)]);
        end
    end
end

% Install required packages if venv exists
if ~isempty(pythonExe) && exist(pythonExe, 'file')
    % Check if torch is already installed
    [status, ~] = system([pythonExe ' -c "import torch" 2>&1']);
    if status ~= 0
        disp('  Installing Python packages (torch, numpy, scipy)...');
        pipCmd = [pythonExe ' -m pip install --quiet --upgrade pip 2>&1'];
        system(pipCmd);
        pipCmd = [pythonExe ' -m pip install --quiet torch numpy scipy 2>&1'];
        [status, result] = system(pipCmd);
        if status == 0
            disp('  Python packages installed successfully.');
        else
            disp(['  Warning: Package installation issue: ' strtrim(result)]);
        end
    else
        disp('  Python packages already installed.');
    end

    % Verify installation
    [~, pyVer] = system([pythonExe ' --version 2>&1']);
    disp(['  Python version: ' strtrim(pyVer)]);

    % Create symlink if venv is not at expected location
    if ~strcmp(venvPath, expectedVenvPath) && ~exist(expectedVenvPath, 'dir')
        disp(['  Creating symlink: ' expectedVenvPath ' -> ' venvPath]);
        % Remove any existing broken symlink
        system(['rm -f "' expectedVenvPath '" 2>&1']);
        [status, result] = system(['ln -sf "' venvPath '" "' expectedVenvPath '" 2>&1']);
        if status ~= 0
            disp(['    Symlink failed: ' strtrim(result)]);
            % Alternative: copy the venv (slower but works)
            disp('    Trying to copy venv instead...');
            system(['cp -r "' venvPath '" "' expectedVenvPath '" 2>&1']);
        end
    end
else
    disp('  Warning: Python venv not available. ProbVer CP-Star may fail.');
end

disp(' ');

% Suppress all warnings to clean up output during test execution
warning('off', 'all');
disp('Warnings suppressed for cleaner output.');
disp(' ');

%% Run Tests
testResults = struct();
testResults.FairNNV = 'not_run';
testResults.ProbVer = 'not_run';
testResults.ModelStar = 'not_run';
testResults.VideoStar = 'not_run';
testResults.GNNV = 'not_run';

% Test 1: ProbVer (run first since it needs Python venv)
disp('========================================');
disp('   TEST 1: ProbVer (Probabilistic Verification)');
disp('========================================');
try
    run_probver();
    testResults.ProbVer = 'passed';
    disp('ProbVer: PASSED');
catch ME
    testResults.ProbVer = 'failed';
    disp(['ProbVer: FAILED - ' ME.message]);
end
disp(' ');

% Test 2: VideoStar
disp('========================================');
disp('   TEST 2: VideoStar');
disp('========================================');
try
    run_videostar();
    testResults.VideoStar = 'passed';
    disp('VideoStar: PASSED');
catch ME
    testResults.VideoStar = 'failed';
    disp(['VideoStar: FAILED - ' ME.message]);
end
disp(' ');

% Test 3: ModelStar
disp('========================================');
disp('   TEST 3: ModelStar (Weight Perturbations)');
disp('========================================');
try
    run_modelstar();
    testResults.ModelStar = 'passed';
    disp('ModelStar: PASSED');
catch ME
    testResults.ModelStar = 'failed';
    disp(['ModelStar: FAILED - ' ME.message]);
end
disp(' ');

% Test 4: GNNV (Graph Neural Network Verification)
disp('========================================');
disp('   TEST 4: GNNV (Graph Neural Network Verification)');
disp('========================================');
try
    run_gnnv();
    testResults.GNNV = 'passed';
    disp('GNNV: PASSED');
catch ME
    testResults.GNNV = 'failed';
    disp(['GNNV: FAILED - ' ME.message]);
end
disp(' ');

% Test 5: FairNNV
disp('========================================');
disp('   TEST 5: FairNNV');
disp('========================================');
try
    run_fairnnv();
    testResults.FairNNV = 'passed';
    disp('FairNNV: PASSED');
catch ME
    testResults.FairNNV = 'failed';
    disp(['FairNNV: FAILED - ' ME.message]);
end
disp(' ');

%% Summary
disp('========================================');
disp('   TEST SUMMARY');
disp('========================================');
disp(['FairNNV:   ' testResults.FairNNV]);
disp(['ProbVer:   ' testResults.ProbVer]);
disp(['ModelStar: ' testResults.ModelStar]);
disp(['VideoStar: ' testResults.VideoStar]);
disp(['GNNV:      ' testResults.GNNV]);
disp(' ');
disp('Results saved to /results/');
disp('========================================');

% Save summary
summaryFile = fullfile(resultsDir, 'test_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'NNV 3.0 Test Suite Summary\n');
fprintf(fid, '==========================\n');
fprintf(fid, 'Date: %s\n\n', datestr(now));
fprintf(fid, 'FairNNV:   %s\n', testResults.FairNNV);
fprintf(fid, 'ProbVer:   %s\n', testResults.ProbVer);
fprintf(fid, 'ModelStar: %s\n', testResults.ModelStar);
fprintf(fid, 'VideoStar: %s\n', testResults.VideoStar);
fprintf(fid, 'GNNV:      %s\n', testResults.GNNV);
fclose(fid);
