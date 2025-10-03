% Test ConcatenationLayer functionality
% To run: results = runtests('test_ConcatenationLayer')

%% Test 1: ConcatenationLayer constructor - with name and dim
L = ConcatenationLayer('concat1', 1);
assert(strcmp(L.Name, 'concat1'));
assert(L.Dim == 1);

%% Test 2: ConcatenationLayer constructor - with all parameters
L = ConcatenationLayer('concat2', 2, 1, {'in1', 'in2'}, {'out'}, 2);
assert(strcmp(L.Name, 'concat2'));
assert(L.Dim == 2);
assert(L.NumInputs == 2);

%% Test 3: ConcatenationLayer evaluate - concatenate along dim 1 (rows)
L = ConcatenationLayer('concat_dim1', 1);

% Create two matrices to concatenate
input1 = [1 2 3; 4 5 6];
input2 = [7 8 9; 10 11 12];

inputs = {input1, input2};
output = L.evaluate(inputs);

% Should concatenate vertically (along rows)
expected = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
assert(isequal(output, expected), 'Concatenation along dim 1 failed');

%% Test 4: ConcatenationLayer evaluate - concatenate along dim 2 (columns)
L = ConcatenationLayer('concat_dim2', 2);

% Create two matrices to concatenate
input1 = [1 2; 3 4];
input2 = [5 6; 7 8];

inputs = {input1, input2};
output = L.evaluate(inputs);

% Should concatenate horizontally (along columns)
expected = [1 2 5 6; 3 4 7 8];
assert(isequal(output, expected), 'Concatenation along dim 2 failed');

%% Test 5: ConcatenationLayer evaluate - three inputs
L = ConcatenationLayer('concat_three', 1);

% Create three matrices
input1 = [1; 2];
input2 = [3; 4];
input3 = [5; 6];

inputs = {input1, input2, input3};
output = L.evaluate(inputs);

% Should concatenate all three vertically
expected = [1; 2; 3; 4; 5; 6];
assert(isequal(output, expected), 'Concatenation of three inputs failed');

%% Test 6: ConcatenationLayer evaluate - 3D arrays along dim 3
L = ConcatenationLayer('concat_3d', 3);

% Create two 3D arrays
input1(:,:,1) = [1 2; 3 4];
input1(:,:,2) = [5 6; 7 8];

input2(:,:,1) = [9 10; 11 12];

inputs = {input1, input2};
output = L.evaluate(inputs);

% Should concatenate along 3rd dimension
assert(size(output, 1) == 2, 'Dim 1 should be unchanged');
assert(size(output, 2) == 2, 'Dim 2 should be unchanged');
assert(size(output, 3) == 3, 'Dim 3 should be sum of inputs');

%% Test 8: ConcatenationLayer reach with ImageStar inputs
L = ConcatenationLayer('concat_imagestar', 3);

% Create two ImageStar inputs
IM1(:,:,1) = [1 1; 0 1];
LB1(:,:,1) = [-0.1 -0.1; 0 0];
UB1(:,:,1) = [0.1 0.1; 0 0];
image_star1 = ImageStar(IM1, LB1, UB1);

IM2(:,:,1) = [0 1; 1 0];
LB2(:,:,1) = [-0.1 -0.1; 0 0];
UB2(:,:,1) = [0.1 0.1; 0 0];
image_star2 = ImageStar(IM2, LB2, UB2);

inputs = {image_star1, image_star2};
output = L.reach(inputs, 'approx-star');

assert(isa(output, 'ImageStar'), 'reach should return ImageStar');

%% Test 9: ConcatenationLayer toGPU
L = ConcatenationLayer('concat', 1);
L_gpu = L.toGPU();
assert(isa(L_gpu, 'ConcatenationLayer'));

%% Test 10: ConcatenationLayer changeParamsPrecision
L = ConcatenationLayer('concat', 1);
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'ConcatenationLayer'));
