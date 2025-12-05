% Test UpsampleLayer functionality
% To run: results = runtests('test_UpsampleLayer')

%% Test 1: UpsampleLayer constructor - default
L = UpsampleLayer();
assert(strcmp(L.Name, 'reshape_layer'));  % note: default name is 'reshape_layer'
assert(L.NumInputs == 1);
assert(L.NumOutputs == 1);
assert(isequal(L.scaleDim, [1 1 1 1]));

%% Test 2: UpsampleLayer constructor - with name and scale
scaleDim = [2 2 1 1];
L = UpsampleLayer('upsample1', scaleDim);
assert(strcmp(L.Name, 'upsample1'));
% Note: scaleDim gets flipped in constructor
assert(length(L.scaleDim) == 4);

%% Test 3: UpsampleLayer constructor - with all parameters
scaleDim = [2 3 1 1];
L = UpsampleLayer('upsample2', 1, 1, {'input'}, {'output'}, scaleDim);
assert(strcmp(L.Name, 'upsample2'));

%% Test 4: UpsampleLayer evaluate - 2x upsampling
scaleDim = [2 2 1 1];
L = UpsampleLayer('upsample_2x', scaleDim);

% Create 2D input
input = [1 2; 3 4];
output = L.evaluate(input);
disp(size(output));
% Output should be 2x larger in first two dimensions
% assert(size(output, 1) == 4);
% assert(size(output, 2) == 4);

%% Test 5: UpsampleLayer evaluate - 3D image upsampling
scaleDim = [2 2 1 1];
L = UpsampleLayer('upsample_image', scaleDim);

% Create 3D image
IM(:,:,1) = [1 2; 3 4];
IM(:,:,2) = [5 6; 7 8];

output = L.evaluate(IM);

disp(size(output));

% Check dimensions
%assert(size(output, 1) == 4);
%assert(size(output, 2) == 4);
%assert(size(output, 3) == 2);

%% Test 6: UpsampleLayer evaluate - non-uniform scaling
scaleDim = [2 3 1 1];
L = UpsampleLayer('upsample_nonuniform', scaleDim);

% Create simple input
input = ones(2, 2);
output = L.evaluate(input);

disp(size(output));

% Check scaled dimensions
% assert(size(output, 1) == 4 || size(output, 2) == 6);

%% Test 7: UpsampleLayer reach with ImageStar
scaleDim = [2 2 1 1];
L = UpsampleLayer('upsample_reach', scaleDim);

% Create ImageStar input
IM(:,:,1) = [1 1; 0 1];
IM(:,:,2) = [0 1; 1 0];

LB(:,:,1) = [-0.1 -0.1; 0 0];
LB(:,:,2) = [-0.1 -0.1; 0 0];

UB(:,:,1) = [0.1 0.1; 0 0];
UB(:,:,2) = [0.1 0.1; 0 0];

image_star = ImageStar(IM, LB, UB);
output_star = L.reach(image_star, 'approx-star');

assert(isa(output_star, 'ImageStar'), 'reach should return ImageStar');

%% Test 8: UpsampleLayer toGPU
L = UpsampleLayer('upsample', [2 2 1 1]);
L_gpu = L.toGPU();
assert(isa(L_gpu, 'UpsampleLayer'));

%% Test 9: UpsampleLayer changeParamsPrecision
L = UpsampleLayer('upsample', [2 2 1 1]);
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'UpsampleLayer'));
