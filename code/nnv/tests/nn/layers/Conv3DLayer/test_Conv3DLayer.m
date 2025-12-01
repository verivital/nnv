% Test Conv3DLayer
% To run: results = runtests('test_Conv3DLayer')

%% Test 1: Conv3DLayer Constructor
% Test basic constructor
% Conv3D weights: [H, W, D, NumChannelsIn, NumFilters]
W = randn(3, 3, 3, 1, 2); % 3x3x3 filter, 1 input channel, 2 filters
b = randn(1, 1, 1, 2); % 2 biases

L = Conv3DLayer(W, b);

assert(isequal(size(L.Weights), [3, 3, 3, 1, 2]));
assert(isequal(size(L.Bias), [1, 1, 1, 2]));

%% Test 2: Conv3DLayer with Parameters
% Test with stride and padding
W = randn(2, 2, 2, 2, 3); % 2x2x2 filter, 2 input channels, 3 filters
b = zeros(1, 1, 1, 3);

name = 'test_conv3d';
padding = [1 1 1; 1 1 1]; % padding on all sides
stride = [1 1 1];
dilation = [1 1 1];

L = Conv3DLayer(name, W, b, padding, stride, dilation);

assert(strcmp(L.Name, name));
assert(isequal(L.PaddingSize, padding));
assert(isequal(L.Stride, stride));
assert(isequal(L.DilationFactor, dilation));

%% Test 3: Conv3DLayer Reach with VolumeStar
% Test reachability analysis
W = ones(2, 2, 2, 1, 1); % Simple filter
b = zeros(1, 1, 1, 1);

L = Conv3DLayer(W, b);

% Create small VolumeStar input
vol_size = [4, 4, 4, 1]; % Small 4x4x4 volume, 1 channel
vol_center = rand(vol_size);
disturbance = 0.1 * ones(vol_size);

inputVolume = VolumeStar(vol_center, -disturbance, disturbance);
outputVolume = L.reach_star_single_input(inputVolume);

assert(isa(outputVolume, 'VolumeStar'));

%% Test 4: Conv3DLayer with Multiple Channels
% Test with multiple input and output channels
numInputChannels = 2;
numFilters = 3;

W = randn(2, 2, 2, numInputChannels, numFilters);
b = randn(1, 1, 1, numFilters);

L = Conv3DLayer(W, b);

vol_size = [4, 4, 4, numInputChannels];
vol_center = rand(vol_size);
disturbance = 0.05 * ones(vol_size);

inputVolume = VolumeStar(vol_center, -disturbance, disturbance);
outputVolume = L.reach_star_single_input(inputVolume);

assert(isa(outputVolume, 'VolumeStar'));

%% Test 5: Conv3DLayer with Stride
% Test with non-unit stride
W = ones(3, 3, 3, 1, 1);
b = zeros(1, 1, 1, 1);

name = 'conv3d_stride';
padding = [0 0 0; 0 0 0];
stride = [2 2 2]; % Stride of 2 in all dimensions
dilation = [1 1 1];

L = Conv3DLayer(name, W, b, padding, stride, dilation);

vol_size = [6, 6, 6, 1];
vol_center = rand(vol_size);
disturbance = 0.1 * ones(vol_size);

inputVolume = VolumeStar(vol_center, -disturbance, disturbance);
outputVolume = L.reach_star_single_input(inputVolume);

assert(isa(outputVolume, 'VolumeStar'));
