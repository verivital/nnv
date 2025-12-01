% test_matlab2nnv
% Unit tests for MATLAB to NNV network conversion
% Tests matlab2nnv function with various network types
% To run: results = runtests('test_matlab2nnv')

%% Test 1: Convert MNIST FC network
rng(42);

% Get path to MNIST model
mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];

% Load network
mnist_model = load(mnist_path);
assert(isfield(mnist_model, 'net'), 'MNIST model should have net field');

% Convert to NNV
net = matlab2nnv(mnist_model.net);

assert(~isempty(net), 'Converted network should not be empty');

%% Test 2: Converted network can evaluate
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);
net = matlab2nnv(mnist_model.net);

% Create a test image (28x28 grayscale)
test_img = single(128 * ones(28, 28));

% Evaluate
output = net.evaluate(test_img);

assert(~isempty(output), 'Output should not be empty');
assert(all(isfinite(output)), 'Output should be finite');

%% Test 3: MNIST output has correct dimensions
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);
net = matlab2nnv(mnist_model.net);

test_img = single(128 * ones(28, 28));
output = net.evaluate(test_img);

% MNIST has 10 classes
assert(length(output) == 10, 'MNIST output should have 10 elements');

%% Test 4: Convert ONNX network via importNetworkFromONNX
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
networks = dir([onnx_path, '*.onnx']);

assert(~isempty(networks), 'Should find ONNX networks');

% Import first network
net_path = [networks(1).folder, filesep, networks(1).name];
matlab_net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');

% Convert to NNV
nnv_net = matlab2nnv(matlab_net);

assert(~isempty(nnv_net), 'Converted ONNX network should not be empty');

%% Test 5: ACAS Xu network can evaluate
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
networks = dir([onnx_path, '*.onnx']);

net_path = [networks(1).folder, filesep, networks(1).name];
matlab_net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');
nnv_net = matlab2nnv(matlab_net);

% ACAS Xu has 5 inputs
test_input = [0.5; 0.5; 0.5; 0.5; 0.5];
output = nnv_net.evaluate(test_input);

assert(length(output) == 5, 'ACAS Xu output should have 5 elements');
assert(all(isfinite(output)), 'Output should be finite');

%% Test 6: Network structure preserved
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);
net = matlab2nnv(mnist_model.net);

% Check that network has layers
assert(net.numLayers > 0, 'Network should have layers');

%% Test 7: Multiple ONNX networks convert
rng(42);

acas_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'NNV2.0', filesep, ...
    'Submission', filesep, 'CAV2023', filesep, 'NNV_vs_MATLAB', filesep, 'acas', filesep];
onnx_path = [acas_path, 'onnx', filesep];
networks = dir([onnx_path, '*.onnx']);

% Convert first 3 networks
for i = 1:min(3, length(networks))
    net_path = [networks(i).folder, filesep, networks(i).name];
    matlab_net = importNetworkFromONNX(net_path, 'InputDataFormats', 'BCSS');
    nnv_net = matlab2nnv(matlab_net);

    assert(~isempty(nnv_net), sprintf('Network %d should convert', i));

    % Test evaluation
    test_input = rand(5, 1);
    output = nnv_net.evaluate(test_input);
    assert(length(output) == 5, sprintf('Network %d should have 5 outputs', i));
end

%% Test 8: Converted network supports reachability
rng(42);

mnist_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'MNIST', filesep, 'mnist_model_fc.mat'];
mnist_model = load(mnist_path);
net = matlab2nnv(mnist_model.net);

% Create a small perturbation ImageStar
img = single(128 * ones(28, 28));
disturbance = 0.5;
lb_clip = max(img - disturbance, 0);
ub_clip = min(img + disturbance, 255);
IS = ImageStar(lb_clip, ub_clip);

% Compute reachable set
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
R = net.reach(IS, reachOptions);

assert(~isempty(R), 'Reachable set should be computed');

