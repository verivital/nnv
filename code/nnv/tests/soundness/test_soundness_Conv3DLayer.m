% test_soundness_Conv3DLayer
% Tests for 3D Convolution layer
% To run: results = runtests('test_soundness_Conv3DLayer')

%% Test 1: Conv3DLayer constructor
rng(42);

% Weights: [H x W x D x NumChannels x NumFilters]
filterH = 3;
filterW = 3;
filterD = 3;
numChannels = 1;
numFilters = 2;

weights = randn(filterH, filterW, filterD, numChannels, numFilters) * 0.1;
bias = randn(1, 1, 1, numFilters) * 0.1;

% Padding must be [2 x 3] matrix: [top bottom; left right; front back] or similar
padding = [0 0 0; 0 0 0];  % 2x3 matrix
stride = [1 1 1];  % 1x3 vector
dilation = [1 1 1];  % 1x3 vector

L = Conv3DLayer('conv3d_test', weights, bias, padding, stride, dilation);

assert(strcmp(L.Name, 'conv3d_test'), 'Name should be set correctly');
assert(L.NumFilters == numFilters, 'NumFilters should match');

%% Test 2: Conv3DLayer with single filter
rng(42);

filterH = 2;
filterW = 2;
filterD = 2;
numChannels = 1;
numFilters = 1;

weights = ones(filterH, filterW, filterD, numChannels, numFilters) * 0.125;
bias = zeros(1, 1, 1, numFilters);

padding = [0 0 0; 0 0 0];
stride = [1 1 1];
dilation = [1 1 1];

L2 = Conv3DLayer('conv3d_single', weights, bias, padding, stride, dilation);

assert(L2.NumFilters == 1, 'Should have 1 filter');

%% Test 3: Conv3DLayer with multiple channels and filters
rng(42);

filterH = 2;
filterW = 2;
filterD = 2;
numChannels = 3;
numFilters = 4;

weights = randn(filterH, filterW, filterD, numChannels, numFilters) * 0.1;
bias = randn(1, 1, 1, numFilters) * 0.1;

padding = [0 0 0; 0 0 0];
stride = [1 1 1];
dilation = [1 1 1];

L3 = Conv3DLayer('conv3d_multi', weights, bias, padding, stride, dilation);

assert(L3.NumFilters == numFilters, 'Should have %d filters', numFilters);
assert(L3.NumChannels == numChannels, 'Should have %d channels', numChannels);

%% Test 4: Conv3DLayer with stride
rng(42);

filterH = 2;
filterW = 2;
filterD = 2;
numChannels = 1;
numFilters = 1;

weights = randn(filterH, filterW, filterD, numChannels, numFilters);
bias = zeros(1, 1, 1, numFilters);

padding = [0 0 0; 0 0 0];
stride = [2 2 2];  % Stride of 2 in all dimensions
dilation = [1 1 1];

L4 = Conv3DLayer('conv3d_stride', weights, bias, padding, stride, dilation);

assert(isequal(L4.Stride, stride), 'Stride should be set correctly');

