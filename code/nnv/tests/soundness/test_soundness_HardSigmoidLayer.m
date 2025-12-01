% test_soundness_HardSigmoidLayer
% Tests for Hard Sigmoid activation layer
% HardSigmoid: f(x) = max(0, min(1, 0.2*x + 0.5))
% To run: results = runtests('test_soundness_HardSigmoidLayer')

%% Test 1: HardSigmoidLayer constructor
rng(42);

L = HardSigmoidLayer('hardsig_test');

assert(strcmp(L.Name, 'hardsig_test'), 'Name should be set');

%% Test 2: HardSigmoidLayer evaluate with known values
rng(42);

L2 = HardSigmoidLayer('hardsig_eval');

% Test known values
% HardSigmoid(0) = 0.2*0 + 0.5 = 0.5
input = [0];
output = L2.evaluate(input);
assert(abs(output - 0.5) < 1e-10, 'HardSigmoid(0) should be 0.5');

%% Test 3: HardSigmoidLayer saturation regions
rng(42);

L3 = HardSigmoidLayer('hardsig_sat');

% Very negative input should saturate to 0
% 0.2*x + 0.5 = 0 => x = -2.5
input_neg = [-10];
output_neg = L3.evaluate(input_neg);
assert(abs(output_neg - 0) < 1e-10, 'Very negative should saturate to 0');

% Very positive input should saturate to 1
% 0.2*x + 0.5 = 1 => x = 2.5
input_pos = [10];
output_pos = L3.evaluate(input_pos);
assert(abs(output_pos - 1) < 1e-10, 'Very positive should saturate to 1');

%% Test 4: HardSigmoidLayer 2D input
rng(42);

L4 = HardSigmoidLayer('hardsig_2d');

input = randn(4, 4);
output = L4.evaluate(input);

% All outputs should be in [0, 1]
assert(all(output(:) >= 0), 'All outputs should be >= 0');
assert(all(output(:) <= 1), 'All outputs should be <= 1');

%% Test 5: HardSigmoidLayer 3D input
% NOTE: HardSig.reach not fully implemented - testing evaluate instead
rng(42);

L5 = HardSigmoidLayer('hardsig_3d');

input = randn(3, 3, 2);
output = L5.evaluate(input);

% All outputs should be in [0, 1]
assert(all(output(:) >= 0), 'All outputs should be >= 0');
assert(all(output(:) <= 1), 'All outputs should be <= 1');

