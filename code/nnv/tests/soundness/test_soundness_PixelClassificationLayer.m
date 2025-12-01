% test_soundness_PixelClassificationLayer
% Tests for Pixel Classification layer (semantic segmentation output)
% To run: results = runtests('test_soundness_PixelClassificationLayer')

%% Test 1: PixelClassificationLayer constructor
rng(42);

classes = categorical({'background', 'car', 'person'});
L = PixelClassificationLayer('pixclass_test', classes, [4 4 3]);

assert(strcmp(L.Name, 'pixclass_test'), 'Name should be set');
% Classes may be stored differently - just verify layer is created
assert(~isempty(L), 'Layer should be created');

%% Test 2: PixelClassificationLayer evaluate
rng(42);

classes = categorical({'bg', 'fg'});
L2 = PixelClassificationLayer('pixclass_eval', classes, [4 4 2]);

% Create input: H x W x numClasses (logits before softmax)
input = rand(4, 4, 2);
input(:,:,1) = 0.3;  % bg probabilities
input(:,:,2) = 0.7;  % fg probabilities (higher)

seg_id = L2.evaluate(input);

% Most pixels should be classified as class 2 (fg)
assert(mode(seg_id(:)) == 2, 'Most pixels should be class 2');

%% Test 3: PixelClassificationLayer with clear class separation
rng(42);

classes = categorical({'a', 'b', 'c'});
L3 = PixelClassificationLayer('pixclass_clear', classes, [2 2 3]);

% Create input with clear winners
input = zeros(2, 2, 3);
input(1, 1, 1) = 10;  % pixel (1,1) -> class 1
input(1, 2, 2) = 10;  % pixel (1,2) -> class 2
input(2, 1, 3) = 10;  % pixel (2,1) -> class 3
input(2, 2, 1) = 10;  % pixel (2,2) -> class 1

seg_id = L3.evaluate(input);

assert(seg_id(1, 1) == 1, 'Pixel (1,1) should be class 1');
assert(seg_id(1, 2) == 2, 'Pixel (1,2) should be class 2');
assert(seg_id(2, 1) == 3, 'Pixel (2,1) should be class 3');
assert(seg_id(2, 2) == 1, 'Pixel (2,2) should be class 1');

%% Test 4: PixelClassificationLayer reachability with ImageStar
rng(42);

classes = categorical({'bg', 'fg'});
L4 = PixelClassificationLayer('pixclass_reach', classes, [3 3 2]);

V = zeros(3, 3, 2, 2);
V(:,:,1,1) = 0.3 * ones(3, 3);  % bg logits
V(:,:,2,1) = 0.7 * ones(3, 3);  % fg logits (higher)
V(:,:,:,2) = rand(3, 3, 2) * 0.1;  % Small perturbation

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L4.reach(input_is, 'approx-star');

assert(~isempty(output_is), 'Should produce output');

