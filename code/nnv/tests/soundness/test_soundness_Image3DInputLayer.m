% test_soundness_Image3DInputLayer
% Tests for 3D Image Input layer
% To run: results = runtests('test_soundness_Image3DInputLayer')

%% Test 1: Image3DInputLayer empty constructor
rng(42);

L = Image3DInputLayer();

assert(strcmp(L.Name, 'Image3dInputLayer'), 'Default name should be set');
assert(strcmp(L.Normalization, 'none'), 'Default normalization should be none');

%% Test 2: Image3DInputLayer with input size
rng(42);

inputSize = [32 32 32 3];  % 32x32x32 with 3 channels
L2 = Image3DInputLayer(inputSize);

assert(isequal(L2.InputSize, inputSize), 'InputSize should be set');

%% Test 3: Image3DInputLayer with normalization
rng(42);

inputSize = [16 16 16 1];
L3 = Image3DInputLayer(inputSize, 'zerocenter');

assert(isequal(L3.InputSize, inputSize), 'InputSize should be set');
assert(strcmp(L3.Normalization, 'zerocenter'), 'Normalization should be zerocenter');

%% Test 4: Image3DInputLayer full constructor
rng(42);

L4 = Image3DInputLayer('img3d_full', [8 8 8 1], 'zerocenter', 'auto', 0.5, 0.1, 0, 1);

assert(strcmp(L4.Name, 'img3d_full'), 'Name should be set');
assert(L4.Mean == 0.5, 'Mean should be set');
assert(L4.StandardDeviation == 0.1, 'Std should be set');

%% Test 5: Image3DInputLayer evaluate (no normalization)
rng(42);

L5 = Image3DInputLayer([4 4 4 1]);

input = rand(4, 4, 4, 1);
output = L5.evaluate(input);

assert(isequal(output, input), 'No normalization should pass through');

