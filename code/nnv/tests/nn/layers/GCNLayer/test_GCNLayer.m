% test GCNLayer class

%% 1) Constructor test - 3 arguments (name, W, b)
W = rand(4, 8);  % 4 input features -> 8 output features
b = rand(8, 1);
L = GCNLayer('test_gcn', W, b);

assert(L.InputSize == 4, 'InputSize should be 4');
assert(L.OutputSize == 8, 'OutputSize should be 8');
assert(strcmp(L.Name, 'test_gcn'), 'Name should match');
assert(isequal(L.Weights, W), 'Weights should match');
assert(isequal(L.Bias, b), 'Bias should match');

%% 2) Constructor test - 2 arguments (W, b)
L2 = GCNLayer(W, b);

assert(L2.InputSize == 4, 'InputSize should be 4');
assert(L2.OutputSize == 8, 'OutputSize should be 8');
assert(strcmp(L2.Name, 'gcn_layer'), 'Default name should be gcn_layer');

%% 3) Constructor test - 0 arguments (empty)
L0 = GCNLayer();

assert(L0.InputSize == 0, 'Empty layer InputSize should be 0');
assert(L0.OutputSize == 0, 'Empty layer OutputSize should be 0');
assert(isempty(L0.Weights), 'Empty layer Weights should be empty');

%% 4) Evaluate test
numNodes = 5;
X = rand(numNodes, 4);  % 5 nodes, 4 features each

% Create a simple normalized adjacency matrix
A = rand(numNodes);
A = A + A';  % symmetrize
D = diag(sum(A, 2));
A_norm = D \ A;  % row-normalize (D^-1 * A)

Y = L.evaluate(X, A_norm);

assert(size(Y, 1) == numNodes, 'Output should have same number of nodes');
assert(size(Y, 2) == 8, 'Output should have 8 features');

% Verify computation manually
expected_Y = A_norm * X * W + b';
assert(max(abs(Y - expected_Y), [], 'all') < 1e-10, 'Evaluate should match manual computation');

%% 5) GraphStar reachability test
NF = rand(numNodes, 4);  % center node features
LB = -0.1 * ones(numNodes, 4);  % perturbation lower bound
UB = 0.1 * ones(numNodes, 4);   % perturbation upper bound

GS_in = GraphStar(NF, LB, UB);

GS_out = L.reach(GS_in, A_norm, 'approx-star');

assert(isa(GS_out, 'GraphStar'), 'Output should be GraphStar');
assert(GS_out.numNodes == numNodes, 'Output should have same number of nodes');
assert(GS_out.numFeatures == 8, 'Output should have 8 features');
assert(GS_out.numPred == GS_in.numPred, 'Number of predicates should be preserved');

%% 6) Verify center matches evaluate
center_in = GS_in.V(:, :, 1);  % center of input GraphStar
center_out = GS_out.V(:, :, 1);  % center of output GraphStar
expected_center = L.evaluate(center_in, A_norm);

assert(max(abs(center_out - expected_center), [], 'all') < 1e-10, 'Center should match evaluate');

%% 7) Verify constraints preserved
assert(isequal(GS_out.C, GS_in.C), 'Constraint matrix C should be preserved');
assert(isequal(GS_out.d, GS_in.d), 'Constraint vector d should be preserved');
assert(isequal(GS_out.pred_lb, GS_in.pred_lb), 'pred_lb should be preserved');
assert(isequal(GS_out.pred_ub, GS_in.pred_ub), 'pred_ub should be preserved');

%% 8) Test without explicit method (should default to approx-star)
GS_out2 = L.reach(GS_in, A_norm);

assert(isa(GS_out2, 'GraphStar'), 'Output should be GraphStar with default method');
assert(max(abs(GS_out2.V - GS_out.V), [], 'all') < 1e-10, 'Default method should give same result');

%% 9) Test precision change
L_single = GCNLayer('single_test', W, b);
L_single.changeParamsPrecision('single');

assert(isa(L_single.Weights, 'single'), 'Weights should be single precision');
assert(isa(L_single.Bias, 'single'), 'Bias should be single precision');

L_single.changeParamsPrecision('double');
assert(isa(L_single.Weights, 'double'), 'Weights should be double precision');

disp('All GCNLayer tests passed!');
