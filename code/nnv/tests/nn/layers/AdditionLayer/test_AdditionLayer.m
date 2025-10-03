% Test AdditionLayer functionality
% To run: results = runtests('test_AdditionLayer')

%% Test 1: AdditionLayer constructor with single name argument
try
    L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});
    assert(strcmp(L.Name, 'add_layer'));
    assert(L.NumInputs == 2);
    assert(L.NumOutputs == 1);
    assert(length(L.InputNames) == 2);
    assert(strcmp(L.InputNames{1}, 'in1'));
    assert(strcmp(L.InputNames{2}, 'in2'));
catch e
    error('AdditionLayer constructor failed: %s', e.message);
end

%% Test 2: AdditionLayer evaluate with two inputs
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});

% Create two simple matrices
input1 = [1 2 3; 4 5 6];
input2 = [1 1 1; 2 2 2];
expected = [2 3 4; 6 7 8];

inputs = {input1, input2};
output = L.evaluate(inputs);

assert(isequal(output, expected), 'AdditionLayer evaluate failed for two inputs');

%% Test 3: AdditionLayer evaluate with three inputs
L = AdditionLayer('add_layer', 3, 1, {'in1', 'in2', 'in3'}, {'out'});

% Create three simple matrices
input1 = [1 2; 3 4];
input2 = [5 6; 7 8];
input3 = [9 10; 11 12];
expected = [15 18; 21 24];

inputs = {input1, input2, input3};
output = L.evaluate(inputs);

assert(isequal(output, expected), 'AdditionLayer evaluate failed for three inputs');

%% Test 4: AdditionLayer evaluate with image data
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});

% Create two 3D image arrays
IM1(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1];
IM1(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1];
IM1(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0];

IM2(:,:,1) = [0 1 0 0; 1 0 0 0; 0 1 0 1; 0 0 0 0];
IM2(:,:,2) = [1 0 1 1; 0 1 1 0; 1 0 0 1; 1 1 1 0];
IM2(:,:,3) = [0 0 0 0; 0 0 1 0; 1 0 0 1; 0 1 0 1];

inputs = {IM1, IM2};
output = L.evaluate(inputs);

expected(:,:,1) = IM1(:,:,1) + IM2(:,:,1);
expected(:,:,2) = IM1(:,:,2) + IM2(:,:,2);
expected(:,:,3) = IM1(:,:,3) + IM2(:,:,3);

assert(isequal(output, expected), 'AdditionLayer evaluate failed for image data');

%% Test 5: AdditionLayer reach with ImageStar inputs
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});

% Create first ImageStar
IM1(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1];
IM1(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1];

LB1(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB1(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

UB1(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB1(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

image_star1 = ImageStar(IM1, LB1, UB1);

% Create second ImageStar
IM2(:,:,1) = [0 1 0 0; 1 0 0 0; 0 1 0 1; 0 0 0 0];
IM2(:,:,2) = [1 0 1 1; 0 1 1 0; 1 0 0 1; 1 1 1 0];

LB2(:,:,1) = [-0.05 -0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB2(:,:,2) = [-0.05 -0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

UB2(:,:,1) = [0.05 0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB2(:,:,2) = [0.05 0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

image_star2 = ImageStar(IM2, LB2, UB2);

% Test reach with approx-star
inputs = {image_star1, image_star2};
output_star = L.reach(inputs, 'approx-star');

assert(isa(output_star, 'ImageStar'), 'reach should return ImageStar');

%% Test 6: AdditionLayer reach with ImageZono inputs
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});

% Create first ImageZono
LB1(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB1(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

UB1(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB1(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

image_zono1 = ImageZono(LB1, UB1);

% Create second ImageZono
LB2(:,:,1) = [-0.05 -0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB2(:,:,2) = [-0.05 -0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

UB2(:,:,1) = [0.05 0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB2(:,:,2) = [0.05 0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

image_zono2 = ImageZono(LB2, UB2);

% Test reach with approx-zono
inputs = {image_zono1, image_zono2};
output_zono = L.reach(inputs, 'approx-zono');

assert(isa(output_zono, 'ImageZono'), 'reach should return ImageZono');

%% Test 7: AdditionLayer toGPU
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});
L_gpu = L.toGPU();
assert(isa(L_gpu, 'AdditionLayer'));

%% Test 8: AdditionLayer changeParamsPrecision
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'AdditionLayer'));
