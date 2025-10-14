% Test SaturatingLinearSymmLayer functionality
% To run: results = runtests('test_SaturatingLinearSymmLayer')
% SatLins (symmetric): f(x) = -1 if x <= -1, x if -1 < x < 1, 1 if x >= 1

%% Test 1: SaturatingLinearSymmLayer constructor - default
L = SaturatingLinearSymmLayer();
assert(strcmp(L.Name, 'act_func_layer'));

%% Test 2: SaturatingLinearSymmLayer constructor - with name
L = SaturatingLinearSymmLayer('satlins1');
assert(strcmp(L.Name, 'satlins1'));

%% Test 3: SaturatingLinearSymmLayer evaluate - simple input
L = SaturatingLinearSymmLayer();

% Test with a simple matrix
input = [-2 -1 -0.5 0 0.5; 1 1.5 2 3 4];
output = L.evaluate(input);

% SatLins: -1 for x<=-1, x for -1<x<1, 1 for x>=1
assert(output(1,1) == -1, 'SatLins(-2) should be -1');
assert(output(1,2) == -1, 'SatLins(-1) should be -1');
assert(output(1,3) == -0.5, 'SatLins(-0.5) should be -0.5');
assert(output(1,4) == 0, 'SatLins(0) should be 0');
assert(output(1,5) == 0.5, 'SatLins(0.5) should be 0.5');
assert(output(2,1) == 1, 'SatLins(1) should be 1');

%% Test 4: SaturatingLinearSymmLayer evaluate - 3D image
L = SaturatingLinearSymmLayer();

% Create 3D image with values in range
IM(:,:,1) = [-1 1 0 -1; 0 0 1 -1; 1 0 -1 0; 1 -1 -1 1];
IM(:,:,2) = [0 1 0 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1];
IM(:,:,3) = [1 -1 1 1; 1 -1 0 1; 0 1 -1 0; 1 0 -1 0];

% Scale to smaller range
IM = 0.5 * IM;

output = L.evaluate(IM);

% All outputs should be in [-1, 1]
assert(all(output(:) >= -1) && all(output(:) <= 1), 'SatLins output should be in [-1,1]');

%% Test 5: SaturatingLinearSymmLayer evaluate - output range check
L = SaturatingLinearSymmLayer();

% Create input with values in middle range (-1, 1)
input = [-0.9 -0.5 -0.1 0 0.1 0.5 0.9];
output = L.evaluate(input);

% Values in (-1,1) should pass through unchanged
assert(all(abs(output - input) < 1e-10), 'Values in (-1,1) should be unchanged');

%% Test 6: SaturatingLinearSymmLayer evaluate - saturation at boundaries
L = SaturatingLinearSymmLayer();

% Test values beyond boundaries
input = [-10; -5; -1; 1; 5; 10];
output = L.evaluate(input);

% Check saturation
assert(output(1) == -1, 'Large negative should saturate at -1');
assert(output(2) == -1, 'Negative beyond -1 should saturate');
assert(output(end-1) == 1, 'Positive beyond 1 should saturate');
assert(output(end) == 1, 'Large positive should saturate at 1');

%% Test 7: SaturatingLinearSymmLayer reach with ImageStar - approx-star
L = SaturatingLinearSymmLayer();

% Create Star input 
lb = [0;0];
ub = [1;1];

I_star = Star(lb,ub);

images = L.reach(I_star, 'approx-star');

assert(~isempty(images), 'SaturatingLinearSymmLayer reach approx-star should return result');
assert(isa(images(1), 'ImageStar') || isa(images(1), 'Star'), 'reach should return ImageStar or Star');

%% Test 8: SaturatingLinearSymmLayer reach with ImageZono
L = SaturatingLinearSymmLayer();

LB(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image_zono = ImageZono(LB, UB);
image = L.reach_zono(image_zono);

assert(isa(image, 'ImageZono') || isa(image, 'Zono'), 'reach_zono should return ImageZono or Zono');

%% Test 9: SaturatingLinearSymmLayer toGPU
L = SaturatingLinearSymmLayer();
L_gpu = L.toGPU();
assert(isa(L_gpu, 'SaturatingLinearSymmLayer'));

%% Test 10: SaturatingLinearSymmLayer changeParamsPrecision
L = SaturatingLinearSymmLayer();
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'SaturatingLinearSymmLayer'));
