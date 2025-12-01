% Test TransposedConv2DLayer (also known as Deconvolution)
% To run: results = runtests('test_TransposedConv2DLayer')

%% Test 1: TransposedConv2DLayer Constructor
% Test basic constructor
% TransposedConv2D weights: [H, W, NumFilters, NumChannels]
% Note: MATLAB drops trailing singleton dimensions, so [3,3,2,1] becomes [3,3,2]
W = randn(3, 3, 2, 1); % 3x3 filter, 2 filters, 1 input channel
b = randn(1, 1, 2); % 2 biases (one per filter)

L = TransposedConv2DLayer(W, b);

% MATLAB drops trailing singleton dim, so [3,3,2,1] -> [3,3,2]
assert(isequal(size(L.Weights), [3, 3, 2]));
assert(isequal(size(L.Bias), [1, 1, 2]));

%% Test 2: TransposedConv2DLayer with Parameters
% Test with stride and cropping (5-arg constructor)
W = ones(4, 4, 3, 2); % 4x4 filter, 3 filters, 2 input channels
b = zeros(1, 1, 3);

name = 'test_transposed_conv';
cropping = [0 0 0 0]; % output cropping [1x4]
stride = [2 2]; % stride for upsampling [1x2]

L = TransposedConv2DLayer(name, W, b, cropping, stride);

assert(strcmp(L.Name, name));
assert(isequal(L.Stride, stride));
assert(isequal(L.CroppingSize, cropping));

%% Test 3: TransposedConv2DLayer Reach
% Test reachability analysis (upsampling)
numChannelsIn = 1;
numChannelsOut = 2;
W = ones(3, 3, numChannelsOut, numChannelsIn);
b = zeros(1, 1, numChannelsOut);

L = TransposedConv2DLayer(W, b);

% Create small ImageStar input
IM = rand(4, 4, numChannelsIn);
LB = -0.1 * ones(4, 4, numChannelsIn);
UB = 0.1 * ones(4, 4, numChannelsIn);

inputImage = ImageStar(IM, LB, UB);
outputImage = L.reach_star_single_input(inputImage);

assert(isa(outputImage, 'ImageStar'));

%% Test 4: TransposedConv2DLayer Upsampling with Stride
% Test upsampling behavior with stride > 1 (4-arg constructor)
W = ones(2, 2, 1, 1);
b = zeros(1, 1, 1);

cropping = [0 0 0 0]; % [1x4]
stride = [2 2]; % 2x upsampling [1x2]

L = TransposedConv2DLayer(W, b, cropping, stride);

IM = ones(3, 3, 1);
LB = -0.05 * ones(3, 3, 1);
UB = 0.05 * ones(3, 3, 1);

inputImage = ImageStar(IM, LB, UB);
outputImage = L.reach_star_single_input(inputImage);

assert(isa(outputImage, 'ImageStar'));

%% Test 5: TransposedConv2DLayer Multiple Channels
% Test with multiple input and output channels
numChannelsIn = 3;
numFilters = 5;
W = randn(3, 3, numFilters, numChannelsIn);
b = randn(1, 1, numFilters);

L = TransposedConv2DLayer(W, b);

IM = rand(5, 5, numChannelsIn);
LB = -0.1 * ones(5, 5, numChannelsIn);
UB = 0.1 * ones(5, 5, numChannelsIn);

inputImage = ImageStar(IM, LB, UB);
outputImage = L.reach_star_single_input(inputImage);

assert(isa(outputImage, 'ImageStar'));
