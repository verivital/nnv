% test_soundness_Conv1DLayer
% Tests for 1D Convolution layer
% To run: results = runtests('test_soundness_Conv1DLayer')

%% Test 1: Conv1DLayer constructor
rng(42);

% Weights: [FilterSize x NumChannels x NumFilters] - must be 3D
filterSize = 3;
numChannels = 2;
numFilters = 4;

weights = randn(filterSize, numChannels, numFilters);
bias = randn(1, numFilters);  % Must be [1 x numFilters]
padding = [0 0];  % [left right]
stride = 1;
dilation = 1;

L = Conv1DLayer('conv1d_test', weights, bias, padding, stride, dilation);

assert(strcmp(L.Name, 'conv1d_test'), 'Name should be set correctly');
assert(L.NumFilters == numFilters, 'NumFilters should match');
assert(L.FilterSize == filterSize, 'FilterSize should match');

%% Test 2: Conv1DLayer with two channels two filters
rng(42);

% NOTE: MATLAB collapses trailing singleton dimensions, so [filterSize, numChannels, 1]
% becomes a 2D array. Must use at least 2 filters to maintain 3D weights.
filterSize = 3;
numChannels = 2;
numFilters = 2;  % Must be >= 2 to avoid MATLAB dimension collapse

weights = randn(filterSize, numChannels, numFilters);
bias = randn(1, numFilters);

L2 = Conv1DLayer('conv1d_simple', weights, bias, [0 0], 1, 1);

assert(L2.NumFilters == 2, 'Should have 2 filters');
assert(L2.NumChannels == 2, 'Should have 2 channels');

%% Test 3: Conv1DLayer with multiple filters
rng(42);

filterSize = 5;
numChannels = 3;
numFilters = 8;

weights = randn(filterSize, numChannels, numFilters) * 0.1;
bias = zeros(1, numFilters);

L3 = Conv1DLayer('conv1d_multi', weights, bias, [0 0], 1, 1);

assert(L3.NumFilters == numFilters, 'Should have %d filters', numFilters);
assert(L3.NumChannels == numChannels, 'Should have %d channels', numChannels);

%% Test 4: Conv1DLayer with padding
rng(42);

filterSize = 3;
numChannels = 1;
numFilters = 2;

weights = randn(filterSize, numChannels, numFilters);
bias = zeros(1, numFilters);
padding = [1 1];  % Same padding

L4 = Conv1DLayer('conv1d_pad', weights, bias, padding, 1, 1);

assert(isequal(L4.PaddingSize, padding), 'Padding should be set correctly');

