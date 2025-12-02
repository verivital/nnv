% test_soundness_SaturatingLinearSymmLayer
% Tests for Symmetric Saturating Linear (Satlins) activation layer
% Satlins: f(x) = -1 if x<-1, x if -1<=x<=1, 1 if x>1
% To run: results = runtests('test_soundness_SaturatingLinearSymmLayer')

%% Test 1: SaturatingLinearSymmLayer constructor
rng(42);

L = SaturatingLinearSymmLayer('satlins_test');

assert(strcmp(L.Name, 'satlins_test'), 'Name should be set');

%% Test 2: SaturatingLinearSymmLayer evaluate with known values
rng(42);

L2 = SaturatingLinearSymmLayer('satlins_eval');

% Test the three regions
% x < -1: output = -1
input_neg = [-2];
output_neg = L2.evaluate(input_neg);
assert(abs(output_neg - (-1)) < 1e-10, 'x < -1 should give -1');

% -1 <= x <= 1: output = x
input_mid = [0.5];
output_mid = L2.evaluate(input_mid);
assert(abs(output_mid - 0.5) < 1e-10, 'Middle region should be identity');

% x > 1: output = 1
input_pos = [2];
output_pos = L2.evaluate(input_pos);
assert(abs(output_pos - 1) < 1e-10, 'x > 1 should give 1');

%% Test 3: SaturatingLinearSymmLayer vector input
rng(42);

L3 = SaturatingLinearSymmLayer('satlins_vec');

input = [-3; -1.5; -1; -0.5; 0; 0.5; 1; 1.5; 3];
expected = [-1; -1; -1; -0.5; 0; 0.5; 1; 1; 1];

output = L3.evaluate(input);
assert(max(abs(output - expected)) < 1e-10, 'Vector evaluation failed');

%% Test 4: SaturatingLinearSymmLayer 2D input
rng(42);

L4 = SaturatingLinearSymmLayer('satlins_2d');

input = randn(4, 4) * 2;  % Some values outside [-1, 1]
output = L4.evaluate(input);

% All outputs should be in [-1, 1]
assert(all(output(:) >= -1), 'All outputs should be >= -1');
assert(all(output(:) <= 1), 'All outputs should be <= 1');

%% Test 5: SaturatingLinearSymmLayer 3D input
% NOTE: SatLins.reach requires Star input - testing evaluate with 3D
rng(42);

L5 = SaturatingLinearSymmLayer('satlins_3d');

input = randn(3, 3, 2) * 2;  % Some values outside [-1, 1]
output = L5.evaluate(input);

% All outputs should be in [-1, 1]
assert(all(output(:) >= -1), 'All outputs should be >= -1');
assert(all(output(:) <= 1), 'All outputs should be <= 1');

