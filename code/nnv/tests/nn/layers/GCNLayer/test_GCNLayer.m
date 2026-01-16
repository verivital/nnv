% test_GCNLayer.m - Unit tests for GCNLayer class
%
% Tests: constructor, evaluate, reach + soundness, precision

% Shared setup (before any %% sections)
W = rand(4, 8); b = rand(8, 1);
L = GCNLayer('test_gcn', W, b);
numNodes = 5;
A_norm = rand(numNodes, numNodes);
X = rand(numNodes, 4);
NF = rand(numNodes, 4);
LB = -0.1 * ones(numNodes, 4);
UB = 0.1 * ones(numNodes, 4);
GS_in = GraphStar(NF, LB, UB);

%% 1) Constructor test
assert(L.InputSize == 4, 'InputSize should be 4');
assert(L.OutputSize == 8, 'OutputSize should be 8');
assert(strcmp(L.Name, 'test_gcn'), 'Name should match');
assert(isequal(L.Weights, W), 'Weights should match');
assert(isequal(L.Bias, b), 'Bias should match');

%% 2) Evaluate test
Y = L.evaluate(X, A_norm);
assert(size(Y, 1) == numNodes, 'Output should have same number of nodes');
assert(size(Y, 2) == 8, 'Output should have 8 features');

% Verify computation manually: Y = A_norm * X * W + b'
expected_Y = A_norm * X * W + b';
assert(max(abs(Y - expected_Y), [], 'all') < 1e-10, 'Evaluate should match manual computation');

%% 3) Reach and soundness test
GS_out = L.reach(GS_in, A_norm, 'approx-star');
assert(isa(GS_out, 'GraphStar'), 'Output should be GraphStar');
assert(GS_out.numNodes == numNodes, 'Output should have same number of nodes');
assert(GS_out.numFeatures == 8, 'Output should have 8 features');
assert(GS_out.numPred == GS_in.numPred, 'Number of predicates should be preserved');

% Soundness: center of output should match evaluate on center of input
center_in = GS_in.V(:, :, 1);
center_out = GS_out.V(:, :, 1);
expected_center = L.evaluate(center_in, A_norm);
assert(max(abs(center_out - expected_center), [], 'all') < 1e-10, 'Center should match evaluate');

% Containment: center output should be within bounds
[lb_out, ub_out] = GS_out.getRanges();
Y_center = L.evaluate(GS_in.V(:,:,1), A_norm);
tol = 1e-6;
assert(all(Y_center(:) >= lb_out(:) - tol), 'Center output should be >= lower bound');
assert(all(Y_center(:) <= ub_out(:) + tol), 'Center output should be <= upper bound');

%% 4) Precision change test
L_prec = GCNLayer('prec_test', W, b);
L_prec.changeParamsPrecision('single');
assert(isa(L_prec.Weights, 'single'), 'Weights should be single precision');
assert(isa(L_prec.Bias, 'single'), 'Bias should be single precision');

L_prec.changeParamsPrecision('double');
assert(isa(L_prec.Weights, 'double'), 'Weights should be double precision');

disp('All GCNLayer tests passed!');
