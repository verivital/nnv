% test_regression_verify_safety
% Regression test for verify_safety method
% To run: results = runtests('test_regression_verify_safety')

%% Test 1: Create network for safety verification
rng(42);
W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
assert(net.numLayers == 2, 'Network should have 2 layers');

%% Test 2: verify_safety method exists
rng(42);
W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
assert(ismethod(net, 'verify_safety'), 'verify_safety method should exist');

%% Test 3: HalfSpace creation for unsafe region
rng(42);
G = [1 0];
g = 5;
H = HalfSpace(G, g);
assert(~isempty(H), 'HalfSpace should be created');
assert(all(H.G == G), 'G matrix should be set');
assert(H.g == g, 'g value should be set');

%% Test 4: verify_safety with clear safe case
rng(42);
W1 = 0.1 * eye(3);
b1 = zeros(3, 1);
W2 = 0.1 * eye(2, 3);
b2 = zeros(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-0.1; -0.1; -0.1];
ub = [0.1; 0.1; 0.1];
I = Star(lb, ub);
G = [1 0];
g = -100;
U = HalfSpace(G, g);
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
[result, ~] = net.verify_safety(I, U, reachOptions, 0);
assert(result == 1, 'Should be safe - output far from unsafe region');

%% Test 5: verify_safety returns valid result
rng(42);
W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-1; -1; -1];
ub = [1; 1; 1];
I = Star(lb, ub);
G = [1 0];
g = 0;
U = HalfSpace(G, g);
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
[result, ~] = net.verify_safety(I, U, reachOptions, 0);
assert(ismember(result, [0, 1, 2]), 'Result should be 0, 1, or 2');

%% Test 6: Multiple HalfSpaces for unsafe region
rng(42);
H1 = HalfSpace([1 0], 5);
H2 = HalfSpace([0 1], 5);
U = [H1; H2];
assert(length(U) == 2, 'Should have 2 HalfSpaces');

%% Test 7: verify_safety with exact-star
rng(42);
W1 = 0.1 * eye(2);
b1 = zeros(2, 1);
W2 = 0.1 * eye(2);
b2 = zeros(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-0.05; -0.05];
ub = [0.05; 0.05];
I = Star(lb, ub);
G = [1 0];
g = -10;
U = HalfSpace(G, g);
reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
[result, ~] = net.verify_safety(I, U, reachOptions, 0);
assert(result == 1, 'Should be safe with exact method');

%% Test 8: verify_safety with falsification samples
rng(42);
W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-0.5; -0.5; -0.5];
ub = [0.5; 0.5; 0.5];
I = Star(lb, ub);
G = [1 0];
g = 0;
U = HalfSpace(G, g);
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
n_samples = 10;
[result, ~] = net.verify_safety(I, U, reachOptions, n_samples);
assert(ismember(result, [0, 1, 2]), 'Result should be valid');
