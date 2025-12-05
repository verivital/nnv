% test_regression_nncs_basic
% Regression test for basic NNCS (Neural Network Control System) operations
% Tests loading controllers, creating plant dynamics, and NNCS objects
% Based on: examples/Tutorial/NNCS/ACC/Verification/verify.m
% To run: results = runtests('test_regression_nncs_basic')

%% Test 1: Load ACC controller network
rng(42);

% Get path to controller
controller_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NNCS', filesep, 'ACC', filesep, 'Verification', filesep, 'controller_5_20.mat'];

% Load controller
net = load_NN_from_mat(controller_path);

assert(~isempty(net), 'Controller network should be loaded');

%% Test 2: Controller network can evaluate
rng(42);

controller_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NNCS', filesep, 'ACC', filesep, 'Verification', filesep, 'controller_5_20.mat'];
net = load_NN_from_mat(controller_path);

% ACC controller takes 5 inputs: [v_set, t_gap, v_ego, D_rel, v_rel]
% Typical values
test_input = [30; 1.4; 32; 80; -2];
output = net.evaluate(test_input);

assert(~isempty(output), 'Controller output should not be empty');
assert(all(isfinite(output)), 'Controller output should be finite');

%% Test 3: Controller output is scalar (acceleration command)
rng(42);

controller_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NNCS', filesep, 'ACC', filesep, 'Verification', filesep, 'controller_5_20.mat'];
net = load_NN_from_mat(controller_path);

test_input = [30; 1.4; 32; 80; -2];
output = net.evaluate(test_input);

% ACC controller outputs a single acceleration command
assert(length(output) == 1, 'ACC controller should output single value');

%% Test 4: Create Star input set
rng(42);

% Initial set for ACC verification
% x = [x_ego, v_ego, a_ego, x_lead, v_lead, a_lead]
lb = [90; 32; 0; 10; 30; 0];
ub = [110; 32.2; 0; 11; 30.2; 0];
init_set = Star(lb, ub);

assert(~isempty(init_set), 'Initial Star set should be created');
assert(init_set.dim == 6, 'Initial set should be 6-dimensional');

%% Test 5: Controller reachability with Star input
rng(42);

controller_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NNCS', filesep, 'ACC', filesep, 'Verification', filesep, 'controller_5_20.mat'];
net = load_NN_from_mat(controller_path);

% Create small input set for controller (5 inputs)
lb = [29; 1.3; 31; 70; -3];
ub = [31; 1.5; 33; 90; -1];
I = Star(lb, ub);

% Compute reachable set
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
R = net.reach(I, reachOptions);

assert(~isempty(R), 'Reachable set should be computed');

%% Test 6: Load controller from Tutorial path (known compatible format)
rng(42);

% Use the Tutorial path which has a compatible format
controller_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NNCS', filesep, 'ACC', filesep, 'Verification', filesep, 'controller_5_20.mat'];

assert(exist(controller_path, 'file') == 2, 'Controller file should exist');

net = load_NN_from_mat(controller_path);
assert(~isempty(net), 'Controller should load');
assert(net.numLayers > 0, 'Controller should have layers');

% Test multiple evaluations
test_inputs = [30 28 32; 1.4 1.2 1.6; 32 30 34; 80 60 100; -2 0 -4];
for i = 1:size(test_inputs, 2)
    output = net.evaluate(test_inputs(:, i));
    assert(isfinite(output), sprintf('Evaluation %d should produce finite output', i));
end

%% Test 7: Controller output bounds reasonable
rng(42);

controller_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NNCS', filesep, 'ACC', filesep, 'Verification', filesep, 'controller_5_20.mat'];
net = load_NN_from_mat(controller_path);

% Create input set
lb = [29; 1.3; 31; 70; -3];
ub = [31; 1.5; 33; 90; -1];
I = Star(lb, ub);

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
R = net.reach(I, reachOptions);

% Get output bounds
[out_lb, out_ub] = R.getRanges;

% Bounds should be reasonable (acceleration typically between -5 and 5 m/s^2)
assert(out_lb >= -10, 'Lower bound should be reasonable');
assert(out_ub <= 10, 'Upper bound should be reasonable');

%% Test 8: Soundness - concrete outputs within reachable bounds
rng(42);

controller_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NNCS', filesep, 'ACC', filesep, 'Verification', filesep, 'controller_5_20.mat'];
net = load_NN_from_mat(controller_path);

lb = [29; 1.3; 31; 70; -3];
ub = [31; 1.5; 33; 90; -1];
I = Star(lb, ub);

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
R = net.reach(I, reachOptions);

[out_lb, out_ub] = R.getRanges;

% Test random samples
n_samples = 20;
tol = 1e-5;
for i = 1:n_samples
    x = lb + (ub - lb) .* rand(5, 1);
    y = net.evaluate(x);
    assert(y >= out_lb - tol && y <= out_ub + tol, ...
        sprintf('Sample %d output should be within bounds', i));
end

