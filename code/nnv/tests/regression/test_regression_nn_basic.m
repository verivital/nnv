% test_regression_nn_basic
% Regression test for basic neural network operations
% Tests NN class creation, evaluation, and reachability
% To run: results = runtests('test_regression_nn_basic')

%% Test 1: Create simple feedforward network
rng(42);

% Create a simple 2-layer network: 3->4->2
W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);

L1 = LayerS(W1, b1, 'poslin');  % ReLU
L2 = LayerS(W2, b2, 'purelin'); % Linear

net = NN({L1; L2});

assert(net.numLayers == 2, 'Network should have 2 layers');

%% Test 2: Network evaluation
rng(42);

W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);

L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});

% Test input
x = [1; 2; 3];
y = net.evaluate(x);

assert(length(y) == 2, 'Output should have 2 elements');
assert(all(isfinite(y)), 'Output should be finite');

% Manual computation for verification
h = max(0, W1 * x + b1);  % ReLU
y_expected = W2 * h + b2; % Linear

assert(max(abs(y - y_expected)) < 1e-10, 'Output should match manual computation');

%% Test 3: Star input reachability (exact)
rng(42);

W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);

L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});

% Input set
lb = [0; 0; 0];
ub = [1; 1; 1];
I = Star(lb, ub);

reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
R = net.reach(I, reachOptions);

assert(~isempty(R), 'Reachable set should be computed');

%% Test 4: Star input reachability (approx)
rng(42);

W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);

L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});

lb = [0; 0; 0];
ub = [1; 1; 1];
I = Star(lb, ub);

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
R = net.reach(I, reachOptions);

assert(~isempty(R), 'Approx reachable set should be computed');
assert(length(R) == 1, 'Approx should return single set');

%% Test 5: Zono input reachability
rng(42);

W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);

L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});

% Create Zono input (approx-zono requires Zono, not Star)
lb = [0; 0; 0];
ub = [1; 1; 1];
c = (lb + ub) / 2;  % center
V = diag((ub - lb) / 2);  % generators
I = Zono(c, V);

reachOptions = struct;
reachOptions.reachMethod = 'approx-zono';
R = net.reach(I, reachOptions);

assert(~isempty(R), 'Zono reachable set should be computed');

%% Test 6: Output bounds correctness
rng(42);

W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);

L1 = LayerS(W1, b1, 'poslin');
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});

lb = [0; 0; 0];
ub = [1; 1; 1];
I = Star(lb, ub);

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
R = net.reach(I, reachOptions);

[out_lb, out_ub] = R.getRanges;

% Test many random points
n_samples = 100;
tol = 1e-5;
for i = 1:n_samples
    x = lb + (ub - lb) .* rand(3, 1);
    y = net.evaluate(x);
    assert(all(y >= out_lb - tol) && all(y <= out_ub + tol), ...
        sprintf('Sample %d output should be within bounds', i));
end

%% Test 7: Deeper network
rng(42);

% 3 hidden layers: 5 -> 8 -> 6 -> 4 -> 3
layers = cell(4, 1);
layers{1} = LayerS(randn(8, 5), randn(8, 1), 'poslin');
layers{2} = LayerS(randn(6, 8), randn(6, 1), 'poslin');
layers{3} = LayerS(randn(4, 6), randn(4, 1), 'poslin');
layers{4} = LayerS(randn(3, 4), randn(3, 1), 'purelin');

net = NN(layers);

assert(net.numLayers == 4, 'Network should have 4 layers');

% Test evaluation
x = randn(5, 1);
y = net.evaluate(x);
assert(length(y) == 3, 'Output should have 3 elements');

%% Test 8: Network with tanh activation
rng(42);

W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);

L1 = LayerS(W1, b1, 'tansig');  % Tanh
L2 = LayerS(W2, b2, 'purelin'); % Linear
net = NN({L1; L2});

x = randn(3, 1);
y = net.evaluate(x);

% Manual computation
h = tanh(W1 * x + b1);
y_expected = W2 * h + b2;

assert(max(abs(y - y_expected)) < 1e-10, 'Tanh network output should be correct');

%% Test 9: Network with sigmoid activation
rng(42);

W1 = randn(4, 3);
b1 = randn(4, 1);
W2 = randn(2, 4);
b2 = randn(2, 1);

L1 = LayerS(W1, b1, 'logsig');  % Sigmoid
L2 = LayerS(W2, b2, 'purelin'); % Linear
net = NN({L1; L2});

x = randn(3, 1);
y = net.evaluate(x);

% Manual computation
h = 1 ./ (1 + exp(-(W1 * x + b1)));
y_expected = W2 * h + b2;

assert(max(abs(y - y_expected)) < 1e-10, 'Sigmoid network output should be correct');

%% Test 10: Soundness - approx-star overapproximates exact-star
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

% Compute exact-star
reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
Re = net.reach(I, reachOptions);

% Union of exact sets
exact_lb = inf(2, 1);
exact_ub = -inf(2, 1);
for i = 1:length(Re)
    [lb_i, ub_i] = Re(i).getRanges;
    exact_lb = min(exact_lb, lb_i);
    exact_ub = max(exact_ub, ub_i);
end

% Compute approx-star
reachOptions.reachMethod = 'approx-star';
Ra = net.reach(I, reachOptions);
[approx_lb, approx_ub] = Ra.getRanges;

% Check that approx-star overapproximates exact-star
tol = 1e-5;
assert(all(approx_lb <= exact_lb + tol), ...
    'approx-star lower bound should overapproximate exact');
assert(all(approx_ub >= exact_ub - tol), ...
    'approx-star upper bound should overapproximate exact');

