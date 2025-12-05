% Test SaturatingLinearLayer functionality
% To run: results = runtests('test_SaturatingLinearLayer')
% SatLin: f(x) = 0 if x <= 0, x if 0 < x < 1, 1 if x >= 1

%% Test 1: SaturatingLinearLayer constructor - default
L = SaturatingLinearLayer();
assert(strcmp(L.Name, 'act_func_layer'));

%% Test 2: SaturatingLinearLayer constructor - with name
L = SaturatingLinearLayer('satlin1');
assert(strcmp(L.Name, 'satlin1'));

%% Test 3: SaturatingLinearLayer evaluate - simple input
L = SaturatingLinearLayer();

% Test with simple matrix
input = [-2 -1 0 0.5 1; 1.5 2 3 4 5];
output = L.evaluate(input);

% SatLin: 0 for x<=0, x for 0<x<1, 1 for x>=1
assert(output(1,1) == 0, 'SatLin(-2) should be 0');
assert(output(1,2) == 0, 'SatLin(-1) should be 0');
assert(output(1,3) == 0, 'SatLin(0) should be 0');
assert(output(1,4) == 0.5, 'SatLin(0.5) should be 0.5');
assert(output(1,5) == 1, 'SatLin(1) should be 1');
assert(output(2,1) == 1, 'SatLin(1.5) should be 1');

%% Test 4: SaturatingLinearLayer evaluate - 3D image
L = SaturatingLinearLayer();

% Create 3D image
IM(:,:,1) = [-1 1 0 -1; 0 0 1 -1; 1 0 -1 0; 1 -1 -1 1];
IM(:,:,2) = [0 1 0 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1];
IM(:,:,3) = [1 -1 1 1; 1 -1 0 1; 0 1 -1 0; 1 0 -1 0];

output = L.evaluate(IM);

% All outputs should be in [0, 1]
assert(all(output(:) >= 0) && all(output(:) <= 1), 'SatLin output should be in [0,1]');

%% Test 5: SaturatingLinearLayer evaluate - output range check
L = SaturatingLinearLayer();

% Create input with values in middle range
input = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
output = L.evaluate(input);

% Values in (0,1) should pass through unchanged
assert(all(abs(output - input) < 1e-10), 'Values in (0,1) should be unchanged');

%% Test 6: SaturatingLinearLayer reach with ImageStar - approx-star
L = SaturatingLinearLayer();

% Create ImageStar input
IM(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1];
IM(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1];
IM(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0];

% Scale to [0,1] range
IM = 0.5 * IM;

LB(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image_star = ImageStar(IM, LB, UB);
images = L.reach(image_star, 'approx-star');

assert(~isempty(images), 'SaturatingLinearLayer reach approx-star should return result');
assert(isa(images(1), 'ImageStar'), 'reach should return ImageStar');

%% Test 7: SaturatingLinearLayer reach with ImageZono
L = SaturatingLinearLayer();

LB(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image_zono = ImageZono(LB, UB);
image = L.reach_zono(image_zono);

assert(isa(image, 'ImageZono'), 'reach_zono should return ImageZono');

%% Test 8: SaturatingLinearLayer reach with Star - approx-star
L = SaturatingLinearLayer();

% Create Star input in valid range
lb = [0.2; 0.3];
ub = [0.7; 0.8];
B = Box(lb, ub);
I_star = B.toStar;

S = L.reach_star_single_input(I_star, 'approx-star', 0, [], 'linprog');
assert(~isempty(S), 'SaturatingLinearLayer reach with Star should return result');

%% Test 9: SaturatingLinearLayer toGPU
L = SaturatingLinearLayer();
L_gpu = L.toGPU();
assert(isa(L_gpu, 'SaturatingLinearLayer'));

%% Test 10: SaturatingLinearLayer changeParamsPrecision
L = SaturatingLinearLayer();
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'SaturatingLinearLayer'));
