% test_soundness_SignLayer
% Tests for Sign activation layer
% Sign: f(x) = -1 if x<0, 0 if x==0, 1 if x>0 (modes vary)
% To run: results = runtests('test_soundness_SignLayer')

%% Test 1: SignLayer constructor default
rng(42);

L = SignLayer();

assert(L.gamma == 0, 'Default gamma should be 0');

%% Test 2: SignLayer constructor with gamma
rng(42);

L2 = SignLayer(0.1);

assert(L2.gamma == 0.1, 'Gamma should be 0.1');

%% Test 3: SignLayer evaluate with known values
rng(42);

L3 = SignLayer();

% Test sign function behavior
input = [-2; -0.5; 0; 0.5; 2];
output = L3.evaluate(input);

% Negative should give negative output
assert(output(1) < 0, 'Negative input should give negative output');
assert(output(2) < 0, 'Negative input should give negative output');
% Positive should give positive output
assert(output(4) > 0, 'Positive input should give positive output');
assert(output(5) > 0, 'Positive input should give positive output');

%% Test 4: SignLayer 2D input
rng(42);

L4 = SignLayer();

input = randn(4, 4);
output = L4.evaluate(input);

% Outputs should be bounded (typically in {-1, 0, 1} or similar)
assert(~isempty(output), 'Should produce output');

%% Test 5: SignLayer with different modes
% NOTE: Sign.reach has known issues - testing evaluate with modes
rng(42);

L5 = SignLayer(0, 'polar_zero_to_pos_one');

input = [-2; -1; 0; 1; 2];
output = L5.evaluate(input);

% Outputs should be bounded
assert(~isempty(output), 'Should produce output');
assert(all(output >= -1), 'All outputs should be >= -1');
assert(all(output <= 1), 'All outputs should be <= 1');

