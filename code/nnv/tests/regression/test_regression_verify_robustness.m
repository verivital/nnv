% test_regression_verify_robustness
% Regression test for neural network robustness verification
% To run: results = runtests('test_regression_verify_robustness')

%% Test 1: Create simple network
rng(42);
W1 = randn(4, 3); b1 = randn(4, 1);
W2 = randn(2, 4); b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
assert(net.numLayers == 2, 'Network should have 2 layers');

%% Test 2: Network evaluation works
rng(42);
W1 = randn(4, 3); b1 = randn(4, 1);
W2 = randn(2, 4); b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
x = randn(3, 1);
y = net.evaluate(x);
assert(length(y) == 2, 'Output should have 2 elements');

%% Test 3: Reachability computes
rng(42);
W1 = randn(4, 3); b1 = randn(4, 1);
W2 = randn(2, 4); b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-0.1; -0.1; -0.1]; ub = [0.1; 0.1; 0.1];
I = Star(lb, ub);
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
R = net.reach(I, reachOptions);
assert(~isempty(R), 'Reachable set should be computed');

%% Test 4: verify_safety method exists
rng(42);
W1 = randn(4, 3); b1 = randn(4, 1);
W2 = randn(2, 4); b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
assert(ismethod(net, 'verify_safety'), 'verify_safety method should exist');

%% Test 5: verify_robustness method exists
rng(42);
W1 = randn(4, 3); b1 = randn(4, 1);
W2 = randn(2, 4); b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
assert(ismethod(net, 'verify_robustness'), 'verify_robustness method should exist');

%% Test 6: verify_vnnlib method exists
rng(42);
W1 = randn(4, 3); b1 = randn(4, 1);
W2 = randn(2, 4); b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
assert(ismethod(net, 'verify_vnnlib'), 'verify_vnnlib method should exist');

%% Test 7: checkRobust method exists
rng(42);
W1 = randn(4, 3); b1 = randn(4, 1);
W2 = randn(2, 4); b2 = randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
assert(ismethod(net, 'checkRobust'), 'checkRobust method should exist');

%% Test 8: Exact-star reachability computes
rng(42);
W1 = 0.5 * randn(3, 2); b1 = 0.1 * randn(3, 1);
W2 = 0.5 * randn(2, 3); b2 = 0.1 * randn(2, 1);
L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-0.05; -0.05]; ub = [0.05; 0.05];
I = Star(lb, ub);
reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
R = net.reach(I, reachOptions);
assert(~isempty(R), 'Exact reachable set should be computed');
