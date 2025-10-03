% Test ReshapeLayer functionality
% To run: results = runtests('test_ReshapeLayer')

%% Test 1: ReshapeLayer constructor - default
L = ReshapeLayer();
assert(strcmp(L.Name, 'reshape_layer'));
assert(L.NumInputs == 1);
assert(L.NumOutputs == 1);

%% Test 2: ReshapeLayer constructor - with name and targetDim
targetDim = [2, 3, 4];
L = ReshapeLayer('reshape1', targetDim);
assert(strcmp(L.Name, 'reshape1'));
assert(isequal(L.targetDim, targetDim));

%% Test 3: ReshapeLayer constructor - with all parameters
targetDim = [4, 4, 1];
L = ReshapeLayer('reshape2', 1, 1, {'input'}, {'output'}, targetDim);
assert(strcmp(L.Name, 'reshape2'));
assert(isequal(L.targetDim, targetDim));

%% Test 4: ReshapeLayer evaluate - simple reshape
targetDim = [2, 3];
L = ReshapeLayer('reshape_test', targetDim);

% Create 1D input
input = [1 2 3 4 5 6];
output = L.evaluate(input);

% Should reshape to [2 3]
assert(all(size(output) == [2 3]), 'ReshapeLayer evaluate failed');
assert(isequal(output, reshape(input, [2 3])), 'Reshape output mismatch');

%% Test 5: ReshapeLayer evaluate - 3D reshape
targetDim = [2, 2, 2];
L = ReshapeLayer('reshape_3d', targetDim);

% Create flat input
input = 1:8;
output = L.evaluate(input);

% Should reshape to [2 2 2]
assert(all(size(output) == [2 2 2]), 'ReshapeLayer 3D evaluate failed');

%% Test 6: ReshapeLayer evaluate - with -1 in targetDim
targetDim = [[], 4];
L = ReshapeLayer('reshape_auto', targetDim);

% Create input
input = [1 2 3 4; 5 6 7 8];
output = L.evaluate(input);

% -1 should be replaced with 1
assert(size(output, 1) == 1 || size(output, 2) == 4, 'ReshapeLayer with -1 failed');

%% Test 7: ReshapeLayer reach with ImageStar
targetDim = [4, 4, 1, 1];
L = ReshapeLayer('reshape_reach', targetDim);

% Create ImageStar input
IM(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1];

LB(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

image_star = ImageStar(IM, LB, UB);
output_star = L.reach(image_star, 'approx-star');

assert(isa(output_star, 'ImageStar'), 'reach should return ImageStar');

%% Test 8: ReshapeLayer toGPU
L = ReshapeLayer('reshape', [2, 3]);
L_gpu = L.toGPU();
assert(isa(L_gpu, 'ReshapeLayer'));

%% Test 9: ReshapeLayer changeParamsPrecision
L = ReshapeLayer('reshape', [2, 3]);
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'ReshapeLayer'));
