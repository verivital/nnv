% test_regression_falsification
% Regression tests for neural network falsification
% Tests NN.falsify method for counterexample finding
% To run: results = runtests('test_regression_falsification')

%% Test 1: falsify method exists
rng(42);
assert(ismethod(NN([]), 'falsify'), 'falsify method should exist');

%% Test 2: Simple network falsification with Star input
rng(42);
W1 = [1 0; 0 1]; b1 = [0; 0];
L1 = LayerS(W1, b1, 'poslin');
W2 = [1 1]; b2 = 0;
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-1; -1]; ub = [1; 1];
I = Star(lb, ub);
G = -1; g = -1.5;
U = HalfSpace(G, g);
n_samples = 100;
counter_inputs = net.falsify(I, U, n_samples);
assert(~isempty(counter_inputs), 'Should find counterexamples');

%% Test 3: Falsification finds valid counterexample
rng(42);
W1 = [1 0; 0 1]; b1 = [0; 0];
L1 = LayerS(W1, b1, 'poslin');
W2 = [1 1]; b2 = 0;
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-1; -1]; ub = [1; 1];
I = Star(lb, ub);
G = -1; g = -1.5;
U = HalfSpace(G, g);
counter_inputs = net.falsify(I, U, 100);
if ~isempty(counter_inputs)
    x = counter_inputs(:,1);
    y = net.evaluate(x);
    assert(U.contains(y), 'Counterexample should reach unsafe region');
    assert(all(x >= lb) && all(x <= ub), 'Counterexample in input set');
end

%% Test 4: Falsification with unreachable unsafe region
rng(42);
W1 = [1 0; 0 1]; b1 = [0; 0];
L1 = LayerS(W1, b1, 'poslin');
W2 = [1 1]; b2 = 0;
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-1; -1]; ub = [1; 1];
I = Star(lb, ub);
G = -1; g = -10;
U = HalfSpace(G, g);
counter_inputs = net.falsify(I, U, 100);
assert(isempty(counter_inputs), 'Should not find counterexamples');

%% Test 5: Falsification with Box input
rng(42);
W1 = [1 0; 0 1]; b1 = [0; 0];
L1 = LayerS(W1, b1, 'poslin');
W2 = [1 1]; b2 = 0;
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [-1; -1]; ub = [1; 1];
I = Box(lb, ub);
G = -1; g = -1.5;
U = HalfSpace(G, g);
counter_inputs = net.falsify(I, U, 100);
assert(~isempty(counter_inputs), 'Should work with Box input');

%% Test 6: Falsification with Zono input directly
rng(42);
W1 = [1 0; 0 1]; b1 = [0; 0];
L1 = LayerS(W1, b1, 'poslin');
W2 = [1 1]; b2 = 0;
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
c = [0; 0]; V = eye(2);
I = Zono(c, V);
G = -1; g = -1.5;
U = HalfSpace(G, g);
counter_inputs = net.falsify(I, U, 100);
assert(~isempty(counter_inputs), 'Should work with Zono input');

%% Test 7: HalfSpace unsafe region specification
rng(42);
G = -1; g = -1.5;
U = HalfSpace(G, g);
assert(isa(U, 'HalfSpace'), 'Unsafe region should be HalfSpace');
assert(U.dim == 1, 'HalfSpace dimension should match output');

%% Test 8: Multiple counterexamples found
rng(42);
W1 = [1 0; 0 1]; b1 = [0; 0];
L1 = LayerS(W1, b1, 'poslin');
W2 = [1 1]; b2 = 0;
L2 = LayerS(W2, b2, 'purelin');
net = NN({L1; L2});
lb = [0; 0]; ub = [2; 2];
I = Star(lb, ub);
G = -1; g = -1;
U = HalfSpace(G, g);
counter_inputs = net.falsify(I, U, 100);
assert(size(counter_inputs, 2) > 1, 'Should find multiple counterexamples');
