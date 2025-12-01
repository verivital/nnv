% test_soundness_TransposedConv1DLayer
% Tests for 1D Transposed Convolution layer
% To run: results = runtests('test_soundness_TransposedConv1DLayer')

%% Test 1: TransposedConv1DLayer constructor (5 args)
rng(42);

% Weights: [FilterSize x NumFilters x NumChannels]
filterSize = 3;
numFilters = 2;
numChannels = 2;

weights = randn(filterSize, numFilters, numChannels);
bias = randn(1, numFilters);
cropping = [0 0];
stride = 1;

L = TransposedConv1DLayer('tconv1d_test', weights, bias, cropping, stride);

assert(strcmp(L.Name, 'tconv1d_test'), 'Name should be set correctly');
assert(L.NumFilters == numFilters, 'NumFilters should match');
assert(L.FilterSize == filterSize, 'FilterSize should match');

%% Test 2: TransposedConv1DLayer with different filter sizes
rng(42);

filterSize = 5;
numFilters = 4;
numChannels = 3;

weights = randn(filterSize, numFilters, numChannels) * 0.1;
bias = zeros(1, numFilters);
cropping = [0 0];
stride = 2;

L2 = TransposedConv1DLayer('tconv1d_stride', weights, bias, cropping, stride);

assert(L2.FilterSize == filterSize, 'FilterSize should be 5');
assert(L2.Stride == stride, 'Stride should be 2');

%% Test 3: TransposedConv1DLayer with cropping
rng(42);

filterSize = 3;
numFilters = 2;
numChannels = 2;

weights = randn(filterSize, numFilters, numChannels);
bias = randn(1, numFilters);
cropping = [1 1];  % Crop 1 from each end
stride = 1;

L3 = TransposedConv1DLayer('tconv1d_crop', weights, bias, cropping, stride);

assert(isequal(L3.CroppingSize, cropping), 'CroppingSize should match');

%% Test 4: TransposedConv1DLayer properties
rng(42);

filterSize = 3;
numFilters = 2;
numChannels = 2;

weights = randn(filterSize, numFilters, numChannels);
bias = randn(1, numFilters);

L4 = TransposedConv1DLayer('tconv1d_props', weights, bias, [0 0], 1);

assert(L4.NumInputs == 1, 'Should have 1 input');
assert(L4.NumOutputs == 1, 'Should have 1 output');
assert(L4.NumChannels == numChannels, 'NumChannels should match');

