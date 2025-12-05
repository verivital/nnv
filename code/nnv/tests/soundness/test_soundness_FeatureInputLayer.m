% test_soundness_FeatureInputLayer
% Tests for Feature Input layer (normalizes input features)
% To run: results = runtests('test_soundness_FeatureInputLayer')

%% Test 1: FeatureInputLayer constructor with InputSize
rng(42);

L = FeatureInputLayer(10);  % 10 features

assert(isequal(L.InputSize, 10), 'InputSize should be 10');

%% Test 2: FeatureInputLayer constructor with normalization
rng(42);

L2 = FeatureInputLayer(5, 'none');

assert(isequal(L2.InputSize, 5), 'InputSize should be 5');
assert(strcmp(L2.Normalization, 'none'), 'Normalization should be none');

%% Test 3: FeatureInputLayer evaluate with no normalization
rng(42);

L3 = FeatureInputLayer(4, 'none');

input = [1; 2; 3; 4];
output = L3.evaluate(input);

% With no normalization, output should equal input
assert(isequal(output, input), 'No normalization should preserve input');

%% Test 4: FeatureInputLayer evaluate with zerocenter normalization
rng(42);

% Create layer with zerocenter normalization
mean_val = [1; 1; 1; 1];
L4 = FeatureInputLayer('feat_zc', 4, 'zerocenter', 'auto', mean_val, [], [], []);

input = [2; 3; 4; 5];
output = L4.evaluate(input);

% Output should be input - mean
expected = input - mean_val;
assert(max(abs(output - expected)) < 1e-10, 'Zerocenter should subtract mean');

%% Test 5: FeatureInputLayer reachability with Star
rng(42);

L5 = FeatureInputLayer(3, 'none');

% Create a simple Star set
c = [1; 2; 3];
V = [c eye(3) * 0.1];
C = [eye(3); -eye(3)];
d = ones(6, 1);

input_star = Star(V, C, d);
output_star = L5.reach(input_star, 'approx-star');

assert(~isempty(output_star), 'Should produce output');
% With no normalization, center should be preserved
assert(max(abs(output_star.V(:,1) - c)) < 1e-10, 'Center should be preserved');

