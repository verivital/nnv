% test_regression_acasxu_vnnlib
% Regression test for ACAS Xu ONNX network with VNN-LIB specification
% Tests ONNX import and VNN-LIB property verification
% Based on: examples/Tutorial/NN/ACAS Xu/verify_onnx_vnnlib.m
% To run: results = runtests('test_regression_acasxu_vnnlib')

%% Test 1: Locate ACAS Xu ONNX networks
rng(42);

% Get path to ACAS Xu data
acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];

% Check ONNX files exist
onnx_path = [acas_path, 'onnx', filesep];
networks = dir([onnx_path, '*.onnx']);

assert(~isempty(networks), 'ACAS Xu ONNX networks should exist');
assert(length(networks) >= 45, 'Should have at least 45 ACAS Xu networks');

%% Test 2: Locate VNN-LIB property files
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];

% Check vnnlib files exist
vnnlib_path = [acas_path, 'vnnlib', filesep];
vnnlibs = dir([vnnlib_path, '*.vnnlib']);

assert(~isempty(vnnlibs), 'VNN-LIB files should exist');
assert(length(vnnlibs) >= 10, 'Should have at least 10 VNN-LIB properties');

%% Test 3: Import ONNX network
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
networks = dir([onnx_path, '*.onnx']);

% Import first network
net_path = [networks(1).folder, filesep, networks(1).name];
net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');

assert(~isempty(net), 'Network should be imported');

%% Test 4: Convert ONNX network to NNV
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
networks = dir([onnx_path, '*.onnx']);

net_path = [networks(1).folder, filesep, networks(1).name];
net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');

% Convert to NNV format
nnv_net = matlab2nnv(net);

assert(~isempty(nnv_net), 'NNV network should be created');

%% Test 5: Load VNN-LIB property
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

% Load property 1
vnnlib_file = [vnnlib_path, 'prop_1.vnnlib'];
property = load_vnnlib(vnnlib_file);

assert(isfield(property, 'lb'), 'Property should have lower bounds');
assert(isfield(property, 'ub'), 'Property should have upper bounds');
assert(isfield(property, 'prop'), 'Property should have output constraints');
assert(length(property.lb) == 5, 'Should have 5 input lower bounds');
assert(length(property.ub) == 5, 'Should have 5 input upper bounds');

%% Test 6: ACAS Xu network evaluation
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
networks = dir([onnx_path, '*.onnx']);

net_path = [networks(1).folder, filesep, networks(1).name];
net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');
nnv_net = matlab2nnv(net);

% Test input (5 features)
test_input = [0.5; 0.5; 0.5; 0.5; 0.5];

% Evaluate
output = nnv_net.evaluate(test_input);

assert(length(output) == 5, 'Output should have 5 elements');
assert(all(isfinite(output)), 'Output should be finite');

%% Test 7: VNN-LIB verification with approx-star
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];
networks = dir([onnx_path, '*.onnx']);

% Use network 1 (1_1)
net_path = [networks(1).folder, filesep, networks(1).name];
net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');
nnv_net = matlab2nnv(net);

% Use property 1
vnnlib_file = [vnnlib_path, 'prop_1.vnnlib'];

% Define reachability options
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

% Verify VNN-LIB property
res = nnv_net.verify_vnnlib(vnnlib_file, reachOptions);

% Result: 0 = unknown, 1 = sat (unsafe), 2 = unsat (safe)
assert(res >= 0 && res <= 2, 'Verification result should be valid');

%% Test 8: HalfSpace property constraint
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];

vnnlib_file = [vnnlib_path, 'prop_1.vnnlib'];
property = load_vnnlib(vnnlib_file);

% Property should contain HalfSpace constraint (G*x <= g)
assert(~isempty(property.prop), 'Property should have output constraints');
assert(isfield(property.prop{1}, 'Hg'), 'Property should have HalfSpace');

hs = property.prop{1}.Hg;
assert(isa(hs, 'HalfSpace'), 'Constraint should be a HalfSpace');

%% Test 9: Verify different property (prop_3)
% NOTE: Direct reach() with Star fails because ONNX network has ImageInputLayer
% which requires ImageStar. Using verify_vnnlib handles this automatically.
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];
networks = dir([onnx_path, '*.onnx']);

net_path = [networks(1).folder, filesep, networks(1).name];
net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');
nnv_net = matlab2nnv(net);

% Verify property 3
vnnlib_file = [vnnlib_path, 'prop_3.vnnlib'];

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
res = nnv_net.verify_vnnlib(vnnlib_file, reachOptions);

assert(res >= 0 && res <= 2, 'Verification result should be valid');

%% Test 10: Verify with different network (network 2)
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];
networks = dir([onnx_path, '*.onnx']);

% Use network 2
net_path = [networks(2).folder, filesep, networks(2).name];
net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');
nnv_net = matlab2nnv(net);

vnnlib_file = [vnnlib_path, 'prop_1.vnnlib'];

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
res = nnv_net.verify_vnnlib(vnnlib_file, reachOptions);

assert(res >= 0 && res <= 2, 'Verification result should be valid');

%% Test 11: Verify property 2
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
vnnlib_path = [acas_path, 'vnnlib', filesep];
networks = dir([onnx_path, '*.onnx']);

net_path = [networks(1).folder, filesep, networks(1).name];
net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');
nnv_net = matlab2nnv(net);

% Verify property 2
vnnlib_file = [vnnlib_path, 'prop_2.vnnlib'];

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
res = nnv_net.verify_vnnlib(vnnlib_file, reachOptions);

assert(res >= 0 && res <= 2, 'Verification result should be valid');

