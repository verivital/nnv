% test_GNN.m - Unit tests for GNN wrapper class
%
% Tests the GNN class functionality including:
% - Constructor variants
% - Forward pass (evaluate)
% - Reachability analysis
% - Graph structure management (setGraph)
%
% Author: Anne Tumlin
% Date: 01/13/2026

%% 1) Empty constructor test
gnn_empty = GNN();
assert(gnn_empty.numLayers == 0, 'Empty GNN should have 0 layers');
assert(isempty(gnn_empty.Layers), 'Empty GNN should have empty Layers');

%% 2) Constructor test - layers only
W1 = rand(4, 8); b1 = rand(8, 1);
W2 = rand(8, 4); b2 = rand(4, 1);
L1 = GCNLayer('gcn1', W1, b1);
L2 = GCNLayer('gcn2', W2, b2);

gnn_layers = GNN({L1, L2});
assert(gnn_layers.numLayers == 2, 'Should have 2 layers');
assert(isempty(gnn_layers.A_norm), 'A_norm should be empty');

%% 3) Constructor test - GCN-only network (layers + A_norm)
numNodes = 5;
A_norm = rand(numNodes, numNodes);
gnn = GNN({L1, L2}, A_norm);

assert(gnn.numLayers == 2, 'Should have 2 layers');
assert(isequal(gnn.A_norm, A_norm), 'A_norm should match');
assert(gnn.InputSize == 4, 'InputSize should be 4');
assert(gnn.OutputSize == 4, 'OutputSize should be 4');

%% 4) Evaluate test - GCN-only
X = rand(numNodes, 4);  % 5 nodes, 4 features
Y = gnn.evaluate(X);

assert(size(Y, 1) == numNodes, 'Output should have same number of nodes');
assert(size(Y, 2) == 4, 'Output should have 4 features');

%% 5) Verify evaluate matches manual layer-by-layer computation
Y_manual = L1.evaluate(X, A_norm);
Y_manual = L2.evaluate(Y_manual, A_norm);
assert(max(abs(Y - Y_manual), [], 'all') < 1e-10, 'GNN.evaluate should match manual computation');

%% 6) Constructor test - Full GNN (layers + A_norm + adj_list + E)
W_node = rand(8, 8); b_node = rand(8, 1);
W_edge = rand(3, 8); b_edge = rand(8, 1);
L_gine = GINELayer('gine', W_node, b_node, W_edge, b_edge);

adj_list = [1 2; 2 3; 3 4; 4 5; 5 1];  % 5 edges forming a cycle
E = rand(5, 3);  % 5 edges, 3 edge features

gnn_full = GNN({L1, L_gine, L2}, A_norm, adj_list, E);

assert(gnn_full.numLayers == 3, 'Should have 3 layers');
assert(isequal(gnn_full.adj_list, adj_list), 'adj_list should match');
assert(isequal(gnn_full.E, E), 'E should match');

%% 7) Constructor test - Full GNN with name
gnn_named = GNN({L1, L2}, A_norm, adj_list, E, [], 'my_gnn');
assert(strcmp(gnn_named.Name, 'my_gnn'), 'Name should match');

%% 8) Evaluate test - Mixed GCN + GINE network
Y_mixed = gnn_full.evaluate(X);

assert(size(Y_mixed, 1) == numNodes, 'Output should have same number of nodes');
assert(size(Y_mixed, 2) == 4, 'Output should have 4 features');

%% 9) Reachability test - GCN-only network
NF = rand(numNodes, 4);
LB = -0.1 * ones(numNodes, 4);
UB = 0.1 * ones(numNodes, 4);
GS_in = GraphStar(NF, LB, UB);

reachOpts = struct('reachMethod', 'approx-star');
GS_out = gnn.reach(GS_in, reachOpts);

assert(isa(GS_out, 'GraphStar'), 'Output should be GraphStar');
assert(GS_out.numNodes == numNodes, 'Output should have same number of nodes');
assert(GS_out.numFeatures == 4, 'Output should have 4 features');

%% 10) Verify center matches evaluate for GCN-only
center_in = GS_in.V(:, :, 1);
center_out = GS_out.V(:, :, 1);
expected = gnn.evaluate(center_in);

assert(max(abs(center_out - expected), [], 'all') < 1e-10, ...
    'Center of output GraphStar should match evaluate()');

%% 11) Verify reachSet and reachTime are populated
assert(length(gnn.reachSet) == gnn.numLayers, 'reachSet should have entry per layer');
assert(length(gnn.reachTime) == gnn.numLayers, 'reachTime should have entry per layer');
assert(all(gnn.reachTime > 0), 'reachTime entries should be positive');

%% 12) Test setGraph - update A_norm only
A_norm_new = rand(numNodes, numNodes);
gnn.setGraph(A_norm_new);

assert(isequal(gnn.A_norm, A_norm_new), 'A_norm should be updated');

%% 13) Test setGraph - weight reuse produces different output
Y_new = gnn.evaluate(X);
assert(~isequal(Y, Y_new), 'Different graph should produce different output');

%% 14) Test setGraph - full update
adj_list_new = [1 3; 2 4; 3 5; 4 1; 5 2];
E_new = rand(5, 3);
gnn_full.setGraph(A_norm_new, adj_list_new, E_new);

assert(isequal(gnn_full.A_norm, A_norm_new), 'A_norm should be updated');
assert(isequal(gnn_full.adj_list, adj_list_new), 'adj_list should be updated');
assert(isequal(gnn_full.E, E_new), 'E should be updated');

%% 15) Test precision change
gnn_prec = GNN({L1, L2}, A_norm);
gnn_prec.changeParamsPrecision('single');

assert(isa(gnn_prec.Layers{1}.Weights, 'single'), 'Weights should be single precision');
assert(isa(gnn_prec.Layers{2}.Weights, 'single'), 'Weights should be single precision');

gnn_prec.changeParamsPrecision('double');
assert(isa(gnn_prec.Layers{1}.Weights, 'double'), 'Weights should be double precision');

%% 16) Test getInfo
info = gnn.getInfo();

assert(info.numLayers == 2, 'numLayers should be 2');
assert(info.hasAdjacency == true, 'hasAdjacency should be true');
assert(strcmp(info.layerTypes{1}, 'GCNLayer'), 'First layer should be GCNLayer');

%% 17) Test with default reachOptions
gnn.setGraph(A_norm);  % Reset to original A_norm
GS_out_default = gnn.reach(GS_in);

assert(isa(GS_out_default, 'GraphStar'), 'Output should be GraphStar with default options');

%% 18) Test output bounds contain samples
num_samples = 10;
[lb_out, ub_out] = GS_out.getRanges();

for s = 1:num_samples
    % Generate random sample from input GraphStar
    alpha = rand(GS_in.numPred, 1) .* (GS_in.pred_ub - GS_in.pred_lb) + GS_in.pred_lb;
    X_sample = GS_in.V(:, :, 1);
    for k = 1:GS_in.numPred
        X_sample = X_sample + alpha(k) * GS_in.V(:, :, k+1);
    end

    % Evaluate at sample
    Y_sample = gnn.evaluate(X_sample);

    % Check sample is within bounds
    tol = 1e-6;
    assert(all(Y_sample(:) >= lb_out(:) - tol), 'Sample should be above lower bound');
    assert(all(Y_sample(:) <= ub_out(:) + tol), 'Sample should be below upper bound');
end

disp('All GNN tests passed!');
