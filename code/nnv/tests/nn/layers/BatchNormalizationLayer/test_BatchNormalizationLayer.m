% Test BatchNormalizationLayer
% To run: results = runtests('test_BatchNormalizationLayer')

%% Test 1: BatchNormalizationLayer Constructor
% Test basic constructor with required parameters
numChannels = 3;
trainedMean = reshape([0.5 0.3 0.7], [1 1 3]); % [1x1xNumChannels]
trainedVariance = reshape([0.1 0.2 0.15], [1 1 3]); % [1x1xNumChannels]

L = BatchNormalizationLayer('NumChannels', numChannels, ...
    'TrainedMean', trainedMean, ...
    'TrainedVariance', trainedVariance);

assert(L.NumChannels == numChannels);
assert(isequal(L.TrainedMean, trainedMean));
assert(isequal(L.TrainedVariance, trainedVariance));
assert(isequal(L.Offset, zeros(1, numChannels))); % Default offset
assert(isequal(L.Scale, ones(1, numChannels))); % Default scale

%% Test 2: BatchNormalizationLayer Constructor with Scale and Offset
% Test constructor with learnable parameters
numChannels = 2;
trainedMean = reshape([1.0 2.0], [1 1 2]);
trainedVariance = reshape([0.5 0.8], [1 1 2]);
scale = reshape([0.9 1.1], [1 1 2]);
offset = reshape([0.1 -0.1], [1 1 2]);

L = BatchNormalizationLayer('NumChannels', numChannels, ...
    'TrainedMean', trainedMean, ...
    'TrainedVariance', trainedVariance, ...
    'Scale', scale, ...
    'Offset', offset);

assert(isequal(L.Scale, scale));
assert(isequal(L.Offset, offset));

%% Test 3: BatchNormalizationLayer Reach with ImageStar
% Test reachability analysis with ImageStar input
numChannels = 2;
trainedMean = reshape([0.5 0.5], [1 1 2]);
trainedVariance = reshape([0.1 0.1], [1 1 2]);
scale = reshape([1.0 1.0], [1 1 2]);
offset = reshape([0.0 0.0], [1 1 2]);

L = BatchNormalizationLayer('NumChannels', numChannels, ...
    'TrainedMean', trainedMean, ...
    'TrainedVariance', trainedVariance, ...
    'Scale', scale, ...
    'Offset', offset);

% Create a simple ImageStar input (2x2 image, 2 channels)
IM = zeros(2, 2, 2);
IM(:, :, 1) = [1.0 2.0; 3.0 4.0];
IM(:, :, 2) = [2.0 3.0; 4.0 5.0];

LB = -0.1 * ones(2, 2, 2);
UB = 0.1 * ones(2, 2, 2);

inputImage = ImageStar(IM, LB, UB);

% Perform reachability analysis
outputImage = L.reach_star_single_input(inputImage);

% Check output is ImageStar
assert(isa(outputImage, 'ImageStar'));

%% Test 4: BatchNormalizationLayer with Single Precision
% Test with single precision inputs
numChannels = 2;
trainedMean = reshape(single([0.5 0.5]), [1 1 2]);
trainedVariance = reshape(single([0.1 0.1]), [1 1 2]);

L = BatchNormalizationLayer('NumChannels', numChannels, ...
    'TrainedMean', trainedMean, ...
    'TrainedVariance', trainedVariance);

% Convert to single precision
L = L.changeParamsPrecision('single');

% Reshape Scale and Offset to 3D (constructor creates 1D by default)
L.Scale = reshape(L.Scale, [1 1 numChannels]);
L.Offset = reshape(L.Offset, [1 1 numChannels]);

% Create single precision input
IM = single(ones(2, 2, 2));
LB = single(-0.1 * ones(2, 2, 2));
UB = single(0.1 * ones(2, 2, 2));

inputImage = ImageStar(IM, LB, UB);
inputImage = inputImage.changeVarsPrecision('single');

% Perform reachability analysis
outputImage = L.reach_star_single_input(inputImage);

% Check output is ImageStar
assert(isa(outputImage, 'ImageStar'));

%% Test 5: BatchNormalizationLayer Multiple Inputs
% Test reach with multiple inputs
numChannels = 1;
trainedMean = reshape(0.0, [1 1 1]);
trainedVariance = reshape(1.0, [1 1 1]);

L = BatchNormalizationLayer('NumChannels', numChannels, ...
    'TrainedMean', trainedMean, ...
    'TrainedVariance', trainedVariance);

% Create two ImageStar inputs
IM1 = ones(3, 3, 1);
IM2 = 2*ones(3, 3, 1);

LB = -0.05 * ones(3, 3, 1);
UB = 0.05 * ones(3, 3, 1);

inputImage1 = ImageStar(IM1, LB, UB);
inputImage2 = ImageStar(IM2, LB, UB);

in_images = [inputImage1 inputImage2];

% Perform reachability with multiple inputs
outputImages = L.reach_star_multipleInputs(in_images, []);

% Check we get two outputs
assert(length(outputImages) == 2);
