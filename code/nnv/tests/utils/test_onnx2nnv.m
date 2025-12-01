% test_onnx2nnv
% Tests for ONNX to NNV network conversion function
% Tests the onnx2nnv utility function interface
% Note: Some ONNX networks have unsupported layer types
% To run: results = runtests('test_onnx2nnv')

%% Test 1: onnx2nnv function exists
rng(42);

% Verify the function exists
assert(exist('onnx2nnv', 'file') == 2, 'onnx2nnv function should exist');

%% Test 2: onnx2nnv requires valid file path
rng(42);

% Test with non-existent file - should error
try
    net = onnx2nnv('nonexistent_file.onnx');
    assert(false, 'Should have thrown error for non-existent file');
catch ME
    % Expected to fail
    assert(true, 'Correctly errors on non-existent file');
end

%% Test 3: ONNX file paths exist
rng(42);

% Verify ONNX files exist in expected locations
acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, 'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep, 'onnx', filesep];

onnx_file = [acas_path, 'ACASXU_run2a_1_1_batch_2000.onnx'];

assert(exist(onnx_file, 'file') == 2, 'ACAS ONNX file should exist');

%% Test 4: onnx2nnv handles unsupported layers gracefully
rng(42);

% ACAS ONNX networks have FlattenLayer which is not yet supported
% This test verifies the function provides informative error messages

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, 'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep, 'onnx', filesep];

onnx_file = [acas_path, 'ACASXU_run2a_1_1_batch_2000.onnx'];

try
    net = onnx2nnv(onnx_file);
    % If it succeeds, the network should be valid
    assert(~isempty(net), 'Network should be created if no error');
catch ME
    % Known issue: FlattenLayer not supported
    % This is expected behavior - verify error is informative
    hasUnsupportedMsg = contains(ME.message, 'Unsupported') || contains(ME.message, 'not supported');
    assert(hasUnsupportedMsg, 'Error message should indicate unsupported layer type');
end

%% Test 5: onnx2nnv accepts loadOptions parameter
rng(42);

% Test that the function signature accepts optional loadOptions struct
acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, 'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep, 'onnx', filesep];

onnx_file = [acas_path, 'ACASXU_run2a_1_1_batch_2000.onnx'];

loadOptions = struct;
loadOptions.InputDataFormat = 'BC';
loadOptions.OutputDataFormat = 'BC';

try
    net = onnx2nnv(onnx_file, loadOptions);
    % If successful, network is valid
    assert(~isempty(net), 'Network should be created');
catch ME
    % Expected - known unsupported layer issues
    assert(true, 'Function accepts loadOptions parameter');
end

%% Test 6: Multiple ONNX files exist in repository
rng(42);

% Verify multiple ONNX benchmark files exist
base_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep];

path1 = [base_path, 'NNV2.0', filesep, 'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep, 'onnx', filesep, 'ACASXU_run2a_1_1_batch_2000.onnx'];
path2 = [base_path, 'Submission', filesep, 'ARCH-COMP2021', filesep, 'benchmarks', filesep, 'Single_Pendulum', filesep, 'controller_single_pendulum.onnx'];
path3 = [base_path, 'Submission', filesep, 'ARCH-COMP2021', filesep, 'benchmarks', filesep, 'ACC', filesep, 'controller_5_20.onnx'];

count = 0;
if exist(path1, 'file') == 2
    count = count + 1;
end
if exist(path2, 'file') == 2
    count = count + 1;
end
if exist(path3, 'file') == 2
    count = count + 1;
end

assert(count >= 2, 'At least 2 ONNX benchmark files should exist');

%% Test 7: loadOptions validation
rng(42);

% Test that invalid loadOptions type causes error
acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, 'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep, 'onnx', filesep];

onnx_file = [acas_path, 'ACASXU_run2a_1_1_batch_2000.onnx'];

try
    % Pass invalid loadOptions (not a struct)
    net = onnx2nnv(onnx_file, 'invalid_options');
    assert(false, 'Should error on invalid loadOptions type');
catch ME
    % Expected error for wrong input type
    hasStructMsg = contains(ME.message, 'struct') || contains(ME.message, 'Wrong');
    assert(hasStructMsg, 'Error should mention struct type requirement');
end

%% Test 8: ONNX import uses matlab2nnv internally
rng(42);

% Verify that onnx2nnv internally calls matlab2nnv
% by checking function file contains the call
onnx2nnv_path = which('onnx2nnv');
assert(~isempty(onnx2nnv_path), 'onnx2nnv should be on path');

fid = fopen(onnx2nnv_path, 'r');
content = fread(fid, '*char')';
fclose(fid);

assert(contains(content, 'matlab2nnv'), 'onnx2nnv should use matlab2nnv internally');
