% Test Conv1DLayer
% To run: results = runtests('test_Conv1DLayer')

%% Test 1: Conv1DLayer Constructor
% Test basic constructor with weights and biases
% Conv1D: [FilterSize, NumChannelsIn, NumFilters]
W = randn(3, 2, 4); % FilterSize=3, 2 input channels, 4 filters
b = randn(1, 4); % 4 biases (must be 2D)

L = Conv1DLayer(W, b);

assert(isequal(size(L.Weights), [3, 2, 4]));
assert(isequal(size(L.Bias), [1, 4]));

%% Test 2: Conv1DLayer with Parameters
% Test with stride, padding, dilation
W = randn(5, 1, 2); % FilterSize=5, 1 input channel, 2 filters
b = zeros(1, 2); % 2D bias

name = 'test_conv1d';
padding = [2 2]; % [1x2] left and right padding
stride = 2; % scalar stride
dilation = 1; % scalar dilation

L = Conv1DLayer(name, W, b, padding, stride, dilation);

assert(strcmp(L.Name, name));
assert(isequal(L.PaddingSize, padding));
assert(isequal(L.Stride, stride));
assert(isequal(L.DilationFactor, dilation));

%% Test 3: Conv1DLayer Reach
% Test reachability analysis
W = ones(3, 1, 2); % Simple filter for testing
b = zeros(1, 2); % 2D bias

L = Conv1DLayer(W, b);

% Create a 1D ImageStar (simulating 1D signal)
% Input: [Length, 1, Channels] format
inputLength = 10;
numChannels = 1;
IM = rand(inputLength, 1, numChannels);
LB = -0.1 * ones(inputLength, 1, numChannels);
UB = 0.1 * ones(inputLength, 1, numChannels);

inputStar = ImageStar(IM, LB, UB);
outputStar = L.reach_star_single_input(inputStar);

assert(isa(outputStar, 'ImageStar'));

%% Test 4: Conv1DLayer Multiple Filters
% Test with multiple input/output channels
numInputChannels = 3;
numFilters = 5;
filterSize = 3;

W = randn(filterSize, numInputChannels, numFilters);
b = randn(1, numFilters); % 2D bias

L = Conv1DLayer(W, b);

inputLength = 16; % Longer input to work with filterSize=3
IM = rand(inputLength, 1, numInputChannels);
LB = -0.05 * ones(inputLength, 1, numInputChannels);
UB = 0.05 * ones(inputLength, 1, numInputChannels);

inputStar = ImageStar(IM, LB, UB);
outputStar = L.reach_star_single_input(inputStar);

assert(isa(outputStar, 'ImageStar'));
