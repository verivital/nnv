% Test HardSigmoidLayer functionality
% To run: results = runtests('test_HardSigmoidLayer')

%% Test 1: HardSigmoidLayer constructor - default
L = HardSigmoidLayer();
assert(strcmp(L.Name, 'act_func_layer'));

%% Test 2: HardSigmoidLayer constructor - with name
L = HardSigmoidLayer('hardsig1');
assert(strcmp(L.Name, 'hardsig1'));

%% Test 3: HardSigmoidLayer evaluate - simple input
L = HardSigmoidLayer();

% Test with a simple matrix
% HardSigmoid: f(x) = max(0, min(1, 0.2*x + 0.5))
input = [-3 -2 -1 0 1; 2 3 4 5 6];
output = L.evaluate(input);

% All outputs should be in [0, 1]
assert(all(output(:) >= 0) && all(output(:) <= 1), 'HardSigmoid output should be in [0,1]');

% Check specific values
% f(-3) = max(0, min(1, -0.6+0.5)) = 0
assert(output(1,1) == 0, 'HardSigmoid(-3) should be 0');

%% Test 4: HardSigmoidLayer evaluate - 3D image
L = HardSigmoidLayer();

% Create 3D image
IM(:,:,1) = [-1 1 0 -1; 0 0 1 -1; 1 0 -1 0; 1 -1 -1 1];
IM(:,:,2) = [0 1 0 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1];
IM(:,:,3) = [1 -1 1 1; 1 -1 0 1; 0 1 -1 0; 1 0 -1 0];

output = L.evaluate(IM);

% All outputs should be in [0, 1]
assert(all(output(:) >= 0) && all(output(:) <= 1), 'HardSigmoid 3D output should be in [0,1]');

%% Test 5: HardSigmoidLayer evaluate - boundary values
L = HardSigmoidLayer();

% Test values at boundaries of the hard sigmoid
input = [-10; -2.5; 0; 2.5; 10];  % outside and at boundaries
output = L.evaluate(input);

% Check saturation at boundaries
assert(output(1) == 0, 'HardSigmoid should saturate at 0 for large negative');
assert(output(end) == 1, 'HardSigmoid should saturate at 1 for large positive');

%% Test 7: HardSigmoidLayer reach with ImageZono
L = HardSigmoidLayer();

LB(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image_zono = ImageZono(LB, UB);
image = L.reach_zono(image_zono);

assert(isa(image, 'ImageZono'), 'reach_zono should return ImageZono');

%% Test 10: HardSigmoidLayer changeParamsPrecision
L = HardSigmoidLayer();
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'HardSigmoidLayer'));
