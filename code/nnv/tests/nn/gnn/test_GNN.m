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

%% 6) SAGEConvLayer reach + soundness
W_node1 = rand(4, 8); W_edge1 = rand(4, 8); b_sage1 = rand(8, 1);
W_node2 = rand(8, 4); W_edge2 = rand(8, 4); b_sage2 = rand(4, 1);
S1 = SAGEConvLayer('sage1', W_node1, W_edge1, b_sage1);
S2 = SAGEConvLayer('sage2', W_node2, W_edge2, b_sage2);

A_binary = double(rand(numNodes, numNodes) > 0.5);
gnn_sage = GNN({S1, S2}, A_binary);

Y_sage = gnn_sage.evaluate(X);
assert(isequal(size(Y_sage), [numNodes, 4]), 'SAGEConv output shape mismatch');

GS_sage = gnn_sage.reach(GS_in, struct('reachMethod', 'approx-star'));
assert(isa(GS_sage, 'GraphStar'), 'SAGEConv reach output should be GraphStar');
center_in_sage = GS_in.V(:, :, 1);
center_out_sage = GS_sage.V(:, :, 1);
assert(max(abs(center_out_sage - gnn_sage.evaluate(center_in_sage)), [], 'all') < 1e-10, ...
    'SAGEConv center mismatch');
[lb_s, ub_s] = GS_sage.getRanges();
Yc_sage = gnn_sage.evaluate(center_in_sage);
assert(all(Yc_sage(:) >= lb_s(:) - tol) && all(Yc_sage(:) <= ub_s(:) + tol), ...
    'SAGEConv center should be within reachable bounds');

%% 7) GINEConvLayer reach + soundness
F_in = 4; F_out = 4; hidden = 8; E_in = 2;
W1g = rand(F_in, hidden); b1g = rand(hidden, 1);
W2g = rand(hidden, F_out); b2g = rand(F_out, 1);
We_g = rand(E_in, F_in); be_g = rand(F_in, 1);
G1 = GINEConvLayer('gineconv1', W1g, b1g, W2g, b2g, We_g, be_g);

% Random small graph: 5 nodes, 8 directed edges
adj_list = [randi(numNodes, 8, 1), randi(numNodes, 8, 1)];
E_feat = rand(size(adj_list, 1), E_in);
edge_w = ones(size(adj_list, 1), 1);
A_norm_g = zeros(numNodes);  % unused for GINE-only stack but required by GNN constructor
gnn_gine = GNN({G1}, A_norm_g, adj_list, E_feat, edge_w);

Y_gine = gnn_gine.evaluate(X);
assert(isequal(size(Y_gine), [numNodes, F_out]), 'GINEConv output shape mismatch');

GS_gine = gnn_gine.reach(GS_in, struct('reachMethod', 'approx-star'));
assert(isa(GS_gine, 'GraphStar'), 'GINEConv reach output should be GraphStar');
center_in_g = GS_in.V(:, :, 1);
center_out_g = GS_gine.V(:, :, 1);
assert(max(abs(center_out_g - gnn_gine.evaluate(center_in_g)), [], 'all') < 1e-8, ...
    'GINEConv center mismatch');
[lb_g, ub_g] = GS_gine.getRanges();
Yc_gine = gnn_gine.evaluate(center_in_g);
assert(all(Yc_gine(:) >= lb_g(:) - tol) && all(Yc_gine(:) <= ub_g(:) + tol), ...
    'GINEConv center should be within reachable bounds');

disp('All GNN tests passed!');
