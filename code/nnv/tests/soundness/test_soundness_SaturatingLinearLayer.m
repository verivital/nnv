% test_soundness_SaturatingLinearLayer
% Tests for Saturating Linear (Satlin) activation layer
% Satlin: f(x) = 0 if x<0, x if 0<=x<=1, 1 if x>1
% To run: results = runtests('test_soundness_SaturatingLinearLayer')

%% Test 1: SaturatingLinearLayer constructor
rng(42);

L = SaturatingLinearLayer('satlin_test');

assert(strcmp(L.Name, 'satlin_test'), 'Name should be set');

%% Test 2: SaturatingLinearLayer evaluate with known values
rng(42);

L2 = SaturatingLinearLayer('satlin_eval');

% Test the three regions
% x < 0: output = 0
input_neg = [-1];
output_neg = L2.evaluate(input_neg);
assert(abs(output_neg - 0) < 1e-10, 'Negative input should give 0');

% 0 <= x <= 1: output = x
input_mid = [0.5];
output_mid = L2.evaluate(input_mid);
assert(abs(output_mid - 0.5) < 1e-10, 'Middle region should be identity');

% x > 1: output = 1
input_pos = [2];
output_pos = L2.evaluate(input_pos);
assert(abs(output_pos - 1) < 1e-10, 'Positive > 1 should give 1');

%% Test 3: SaturatingLinearLayer vector input
rng(42);

L3 = SaturatingLinearLayer('satlin_vec');

input = [-2; -0.5; 0; 0.5; 1; 1.5; 3];
expected = [0; 0; 0; 0.5; 1; 1; 1];

output = L3.evaluate(input);
assert(max(abs(output - expected)) < 1e-10, 'Vector evaluation failed');

%% Test 4: SaturatingLinearLayer 2D input
rng(42);

L4 = SaturatingLinearLayer('satlin_2d');

input = randn(4, 4);
output = L4.evaluate(input);

% All outputs should be in [0, 1]
assert(all(output(:) >= 0), 'All outputs should be >= 0');
assert(all(output(:) <= 1), 'All outputs should be <= 1');

%% Test 5: SaturatingLinearLayer reachability with ImageStar
rng(42);

L5 = SaturatingLinearLayer('satlin_reach');

V = zeros(3, 3, 1, 2);
V(:,:,1,1) = 0.5 * ones(3, 3);  % Center at 0.5
V(:,:,1,2) = rand(3, 3) * 0.2;  % Small perturbation

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L5.reach(input_is, 'approx-star');

assert(~isempty(output_is), 'Should produce output');

