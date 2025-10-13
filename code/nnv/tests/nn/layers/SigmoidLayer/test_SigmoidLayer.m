% Test SigmoidLayer functionality
% To run: results = runtests('test_SigmoidLayer')

%% Test 1: SigmoidLayer constructor - default
L = SigmoidLayer();
assert(strcmp(L.Name, 'act_func_layer'));

%% Test 2: SigmoidLayer constructor - with name
L = SigmoidLayer('sigmoid1');
assert(strcmp(L.Name, 'sigmoid1'));

%% Test 3: SigmoidLayer evaluate - simple input
L = SigmoidLayer();

% Test with simple matrix
input = [-2 -1 0 1 2; -3 -0.5 0.5 3 4];
output = L.evaluate(input);

% Sigmoid should be applied element-wise
expected = 1 ./ (1 + exp(-input));
assert(all(abs(output(:) - expected(:)) < 1e-6), 'SigmoidLayer evaluate failed');

%% Test 4: SigmoidLayer evaluate - 3D image
L = SigmoidLayer();

% Create 3D image
IM(:,:,1) = [-1 1 0 -1; 0 0 1 -1; 1 0 -1 0; 1 -1 -1 1];
IM(:,:,2) = [0 1 0 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1];
IM(:,:,3) = [1 -1 1 1; 1 -1 0 1; 0 1 -1 0; 1 0 -1 0];

output = L.evaluate(IM);

% Sigmoid should be applied element-wise
expected = 1 ./ (1 + exp(-IM));
assert(all(abs(output(:) - expected(:)) < 1e-6), 'SigmoidLayer 3D evaluate failed');

%% Test 5: SigmoidLayer evaluate - output range check
L = SigmoidLayer();

% Create input with wide range
input = [-10 -5 0 5 10];
output = L.evaluate(input);

% All sigmoid outputs should be in (0, 1)
assert(all(output(:) > 0) && all(output(:) < 1), 'Sigmoid output should be in (0,1)');

%% Test 6: SigmoidLayer reach with ImageStar - approx-star
L = SigmoidLayer();

% Create ImageStar input
IM(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1];
IM(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1];
IM(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0];

LB(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image_star = ImageStar(IM, LB, UB);
images = L.reach(image_star, 'approx-star');

assert(~isempty(images), 'SigmoidLayer reach approx-star should return result');
assert(isa(images(1), 'ImageStar'), 'reach should return ImageStar');

%% Test 7: SigmoidLayer reach with ImageZono
L = SigmoidLayer();

LB(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image_zono = ImageZono(LB, UB);
image = L.reach_zono(image_zono);

assert(isa(image, 'ImageZono'), 'reach_zono should return ImageZono');

%% Test 8: SigmoidLayer reach with Star - approx-star
L = SigmoidLayer();

% Create Star input
lb = [-1; -1];
ub = [1; 1];
B = Box(lb, ub);
I_star = B.toStar;

S = L.reach_star_single_input(I_star, 'approx-star', 0, [], 'linprog');
assert(~isempty(S), 'SigmoidLayer reach with Star should return result');

%% Test 9: SigmoidLayer toGPU
L = SigmoidLayer();
L_gpu = L.toGPU();
assert(isa(L_gpu, 'SigmoidLayer'));

%% Test 10: SigmoidLayer changeParamsPrecision
L = SigmoidLayer();
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'SigmoidLayer'));
