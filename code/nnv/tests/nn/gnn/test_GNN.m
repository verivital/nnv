% test_GNN.m - Unit tests for GNN wrapper class
%
% Tests: constructor, evaluate, reach + soundness, setGraph, precision

% Shared setup (before any %% sections)
W1 = rand(4, 8); b1 = rand(8, 1);
W2 = rand(8, 4); b2 = rand(4, 1);
L1 = GCNLayer('gcn1', W1, b1);
L2 = GCNLayer('gcn2', W2, b2);

numNodes = 5;
A_norm = rand(numNodes, numNodes);
X = rand(numNodes, 4);

NF = rand(numNodes, 4);
LB = -0.1 * ones(numNodes, 4);
UB = 0.1 * ones(numNodes, 4);
GS_in = GraphStar(NF, LB, UB);

gnn = GNN({L1, L2}, A_norm);

%% 1) Constructor test
assert(gnn.numLayers == 2, 'Should have 2 layers');
assert(isequal(gnn.A_norm, A_norm), 'A_norm should match');
assert(gnn.InputSize == 4, 'InputSize should be 4');
assert(gnn.OutputSize == 4, 'OutputSize should be 4');

%% 2) Evaluate test
Y = gnn.evaluate(X);
assert(size(Y, 1) == numNodes, 'Output should have same number of nodes');
assert(size(Y, 2) == 4, 'Output should have 4 features');

% Verify matches manual layer-by-layer computation
Y_manual = L1.evaluate(X, A_norm);
Y_manual = L2.evaluate(Y_manual, A_norm);
assert(max(abs(Y - Y_manual), [], 'all') < 1e-10, 'GNN.evaluate should match manual computation');

%% 3) Reach and soundness test
reachOpts = struct('reachMethod', 'approx-star');
GS_out = gnn.reach(GS_in, reachOpts);

assert(isa(GS_out, 'GraphStar'), 'Output should be GraphStar');
assert(GS_out.numNodes == numNodes, 'Output should have same number of nodes');
assert(GS_out.numFeatures == 4, 'Output should have 4 features');

% Soundness: center should match evaluate
center_in = GS_in.V(:, :, 1);
center_out = GS_out.V(:, :, 1);
expected = gnn.evaluate(center_in);
assert(max(abs(center_out - expected), [], 'all') < 1e-10, ...
    'Center of output GraphStar should match evaluate()');

% Containment: center output should be within bounds
[lb_out, ub_out] = GS_out.getRanges();
Y_center = gnn.evaluate(GS_in.V(:,:,1));
tol = 1e-6;
assert(all(Y_center(:) >= lb_out(:) - tol), 'Center output should be >= lower bound');
assert(all(Y_center(:) <= ub_out(:) + tol), 'Center output should be <= upper bound');

% Verify reachSet and reachTime populated
assert(length(gnn.reachSet) == gnn.numLayers, 'reachSet should have entry per layer');
assert(length(gnn.reachTime) == gnn.numLayers, 'reachTime should have entry per layer');
assert(all(gnn.reachTime > 0), 'reachTime entries should be positive');

%% 4) setGraph test
Y_original = gnn.evaluate(X);  % Store original output
A_norm_new = rand(numNodes, numNodes);
gnn.setGraph(A_norm_new);
assert(isequal(gnn.A_norm, A_norm_new), 'A_norm should be updated');

Y_new = gnn.evaluate(X);
assert(~isequal(Y_original, Y_new), 'Different graph should produce different output');

%% 5) Precision change test
gnn_prec = GNN({L1, L2}, A_norm);
gnn_prec.changeParamsPrecision('single');
assert(isa(gnn_prec.Layers{1}.Weights, 'single'), 'Weights should be single precision');

gnn_prec.changeParamsPrecision('double');
assert(isa(gnn_prec.Layers{1}.Weights, 'double'), 'Weights should be double precision');

disp('All GNN tests passed!');
