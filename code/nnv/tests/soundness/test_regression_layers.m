% test_regression_layers
% Regression tests with known inputs and expected outputs
% These tests capture expected behavior to detect regressions
% To run: results = runtests('test_regression_layers')

%% Test 1: ReLU with known values
rng(42);

L = ReluLayer();

% Known input -> expected output
input = [-2; -1; 0; 1; 2];
expected = [0; 0; 0; 1; 2];

output = L.evaluate(input);
assert(isequal(output, expected), 'ReLU regression failed');

%% Test 2: LeakyReLU with known values
rng(42);

gamma = 0.1;
L2 = LeakyReluLayer('leaky_reg', 1, {'in'}, 1, {'out'}, gamma);

input = [-2; -1; 0; 1; 2];
expected = [-0.2; -0.1; 0; 1; 2];

output = L2.evaluate(input);
assert(max(abs(output - expected)) < 1e-10, 'LeakyReLU regression failed');

%% Test 3: Sigmoid with known values
rng(42);

L3 = SigmoidLayer('sig_reg');

input = [0];
expected = [0.5];

output = L3.evaluate(input);
assert(abs(output - expected) < 1e-10, 'Sigmoid(0) should be 0.5');

%% Test 4: Tanh with known values
rng(42);

L4 = TanhLayer('tanh_reg');

input = [0];
expected = [0];

output = L4.evaluate(input);
assert(abs(output - expected) < 1e-10, 'Tanh(0) should be 0');

%% Test 5: FullyConnected with known weights
rng(42);

W = [1 2; 3 4];  % 2x2 matrix
b = [0.5; 1.0];  % 2x1 bias

L5 = FullyConnectedLayer(W, b);

input = [1; 1];
% Expected: W * input + b = [1*1+2*1+0.5; 3*1+4*1+1.0] = [3.5; 8.0]
expected = [3.5; 8.0];

output = L5.evaluate(input);
assert(max(abs(output - expected)) < 1e-10, 'FC regression failed');

%% Test 6: Conv2D with known filter
rng(42);

% Use 2 filters to avoid MATLAB dimension collapse with single filter
W = randn(3, 3, 1, 2) * 0.1;  % 2 filters
b = zeros(1, 1, 2);  % 3D bias: [1, 1, numFilters]

L6 = Conv2DLayer(W, b, [0 0 0 0], [1 1], [1 1]);

input = rand(5, 5);
output = L6.evaluate(input);

% Just verify output is produced
assert(~isempty(output), 'Conv2D should produce output');

%% Test 7: MaxPooling2D with known values
rng(42);

L7 = MaxPooling2DLayer('maxpool_reg', [2 2], [2 2], [0 0 0 0]);

input = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];
% 2x2 pools: max([1,2,5,6])=6, max([3,4,7,8])=8, max([9,10,13,14])=14, max([11,12,15,16])=16
expected = [6 8; 14 16];

output = L7.evaluate(input);
assert(isequal(output, expected), 'MaxPool regression failed');

%% Test 8: AveragePooling2D with known values
rng(42);

L8 = AveragePooling2DLayer('avgpool_reg', [2 2], [2 2], [0 0 0 0]);

input = [2 4 6 8; 2 4 6 8; 2 4 6 8; 2 4 6 8];
% 2x2 pools: mean([2,4,2,4])=3, mean([6,8,6,8])=7, etc.
expected = [3 7; 3 7];

output = L8.evaluate(input);
assert(max(abs(output(:) - expected(:))) < 1e-10, 'AvgPool regression failed');

%% Test 9: FlattenLayer with known shape
rng(42);

L9 = FlattenLayer('flatten_reg');
L9.Type = 'nnet.cnn.layer.FlattenLayer';  % Set required Type

input = rand(2, 3, 4);  % 24 elements
output = L9.evaluate(input);

assert(numel(output) == 24, 'Flatten should preserve element count');

%% Test 10: BatchNormalization with known parameters
rng(42);

% Simple case: mean=0, var=1, scale=1, offset=0 should be identity
L10 = BatchNormalizationLayer('bn_reg', 0, 1, 0, 1, 1);

input = [1 2 3 4 5];
expected = input;  % Identity

output = L10.evaluate(input');
assert(max(abs(output' - expected)) < 1e-5, 'BatchNorm identity regression failed');

