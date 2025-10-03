% Test PixelClassificationLayer functionality
% To run: results = runtests('test_PixelClassificationLayer')

%% Test 1: PixelClassificationLayer constructor - 3 parameters
classes = categorical({'background', 'object'});
outputSize = [32 32 2];
L = PixelClassificationLayer('pixelclass1', classes, outputSize);

assert(strcmp(L.Name, 'pixelclass1'));
assert(isequal(L.OutputSize, outputSize));
assert(L.NumInputs == 1);

%% Test 2: PixelClassificationLayer constructor - verify classes
classes = categorical({'class1', 'class2', 'class3'});
outputSize = [16 16 3];
L = PixelClassificationLayer('pixelclass2', classes, outputSize);

% Should have original classes plus 'unknown' and 'misclass'
assert(length(categories(L.Classes)) == 5, 'Should have 5 classes (3 + unknown + misclass)');

%% Test 3: PixelClassificationLayer getClasses
classes = categorical({'background', 'object'});
outputSize = [32 32 2];
L = PixelClassificationLayer('pixelclass3', classes, outputSize);

% Get classes for specific indices
idxs = [1 2];
[result_id, result_classes] = L.getClasses(idxs);
assert(isa(result_classes, 'categorical'), 'Should return categorical array');

%% Test 4: PixelClassificationLayer classify - simple
classes = categorical({'background', 'foreground'});
outputSize = [2 2 2];
L = PixelClassificationLayer('pixelclass4', classes, outputSize);

% Create simple segmentation map (2x2 with 2 classes)
seg_map(:,:,1) = [0.8 0.2; 0.3 0.7];  % background scores
seg_map(:,:,2) = [0.2 0.8; 0.7 0.3];  % foreground scores

[seg_result, class_names] = L.evaluate(seg_map);

% Check dimensions
assert(size(seg_result, 1) == 2, 'Height should be 2');
assert(size(seg_result, 2) == 2, 'Width should be 2');

%% Test 5: PixelClassificationLayer evaluateSegmentation
classes = categorical({'class1', 'class2'});
outputSize = [4 4 2];
L = PixelClassificationLayer('pixelclass5', classes, outputSize);

% Create ImageStar representing segmentation scores
IM(:,:,1) = ones(4, 4);  % class1 scores
IM(:,:,2) = zeros(4, 4);  % class2 scores
LB = -0.1 * ones(4, 4, 2);
UB = 0.1 * ones(4, 4, 2);

image_star = ImageStar(IM, LB, UB);
seg_result = L.reach(image_star,'approx-star');

% Should classify based on maximum scores
assert(~isempty(seg_result), 'Should return segmentation result');

%% Test 6: PixelClassificationLayer reach
classes = categorical({'background', 'object'});
outputSize = [2 2 2];
L = PixelClassificationLayer('pixelclass6', classes, outputSize);

% Create ImageStar input
IM(:,:,1) = [1 1; 0 1];
IM(:,:,2) = [0 0; 1 0];
LB = -0.1 * ones(2, 2, 2);
UB = 0.1 * ones(2, 2, 2);

image_star = ImageStar(IM, LB, UB);
output = L.reach(image_star, 'approx-star');

% Pixel classification layer should pass through unchanged
assert(isa(output, 'Star'), 'reach should return ImageStar');

%% Test 7: PixelClassificationLayer toGPU
classes = categorical({'class1', 'class2'});
outputSize = [4 4 2];
L = PixelClassificationLayer('pixelclass7', classes, outputSize);

L_gpu = L.toGPU();
assert(isa(L_gpu, 'PixelClassificationLayer'));

%% Test 8: PixelClassificationLayer changeParamsPrecision
classes = categorical({'class1', 'class2'});
outputSize = [4 4 2];
L = PixelClassificationLayer('pixelclass8', classes, outputSize);

L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'PixelClassificationLayer'));
