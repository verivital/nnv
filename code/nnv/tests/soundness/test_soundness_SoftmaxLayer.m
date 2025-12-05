% test_soundness_SoftmaxLayer
% Tests for Softmax layer
% To run: results = runtests('test_soundness_SoftmaxLayer')

%% Test 1: SoftmaxLayer constructor
rng(42);

L = SoftmaxLayer('softmax_test');

assert(strcmp(L.Name, 'softmax_test'), 'Name should be set');

%% Test 2: SoftmaxLayer evaluate
rng(42);

L2 = SoftmaxLayer('softmax_eval');

input = randn(5, 1);  % 5 classes
output = L2.evaluate(input);

% Softmax output should sum to 1
assert(abs(sum(output) - 1) < 1e-10, 'Softmax outputs should sum to 1');

% All outputs should be positive
assert(all(output > 0), 'All softmax outputs should be positive');

%% Test 3: SoftmaxLayer with known values
rng(42);

L3 = SoftmaxLayer('softmax_known');

% Simple test: [0, 0, 0] -> [1/3, 1/3, 1/3]
input = zeros(3, 1);
output = L3.evaluate(input);

expected = ones(3, 1) / 3;
assert(max(abs(output - expected)) < 1e-10, 'Uniform input should give uniform output');

%% Test 4: SoftmaxLayer with extreme values
rng(42);

L4 = SoftmaxLayer('softmax_extreme');

% Large positive value should dominate
input = [10; 0; 0];
output = L4.evaluate(input);

assert(output(1) > 0.99, 'First output should be close to 1');
assert(sum(output) - 1 < 1e-10, 'Sum should still be 1');

%% Test 5: SoftmaxLayer 3D image input
rng(42);

L5 = SoftmaxLayer('softmax_3d');

% Use 3D image format (H x W x C)
input = randn(1, 1, 4);  % 1x1 spatial, 4 channels
output = L5.evaluate(input);

% Output should sum to 1 across channels
assert(abs(sum(output(:)) - 1) < 1e-10, 'Softmax outputs should sum to 1');
assert(all(output(:) > 0), 'All softmax outputs should be positive');

