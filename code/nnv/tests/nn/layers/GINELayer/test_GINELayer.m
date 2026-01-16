% test_GINELayer.m - Unit tests for GINELayer class
%
% Tests: constructor, evaluate, reach + soundness, edge perturbation, precision

% Shared setup (before any %% sections)
W_node = rand(4, 8); b_node = rand(8, 1);
W_edge = rand(3, 8); b_edge = rand(8, 1);
L = GINELayer('test_gine', W_node, b_node, W_edge, b_edge);

numNodes = 4;
numEdges = 5;
X = rand(numNodes, 4);
E = rand(numEdges, 3);
adj_list = [1 2; 1 3; 2 3; 3 4; 4 1];

NF = rand(numNodes, 4);
LB = -0.1 * ones(numNodes, 4);
UB = 0.1 * ones(numNodes, 4);
GS_in = GraphStar(NF, LB, UB);

%% 1) Constructor test
assert(L.InputSize == 4, 'InputSize should be 4');
assert(L.OutputSize == 8, 'OutputSize should be 8');
assert(L.EdgeInputSize == 3, 'EdgeInputSize should be 3');
assert(strcmp(L.Name, 'test_gine'), 'Name should match');
assert(isequal(L.Weights, W_node), 'Node weights should match');
assert(isequal(L.EdgeWeights, W_edge), 'Edge weights should match');
assert(L.Epsilon == 0, 'Default epsilon should be 0');

%% 2) Evaluate test
Y = L.evaluate(X, E, adj_list);
assert(size(Y, 1) == numNodes, 'Output should have same number of nodes');
assert(size(Y, 2) == 8, 'Output should have 8 features');

% Verify computation manually (GNNV architecture)
H = X * W_node;  % Transform nodes first
src_nodes = adj_list(:, 1);
dst_nodes = adj_list(:, 2);
H_src = H(src_nodes, :);
E_trans = E * W_edge + b_edge';
edge_msg = max(0, H_src + E_trans);  % ReLU
agg = zeros(numNodes, 8);
for e = 1:numEdges
    agg(dst_nodes(e), :) = agg(dst_nodes(e), :) + edge_msg(e, :);
end
expected_Y = H + agg + b_node';
assert(max(abs(Y - expected_Y), [], 'all') < 1e-10, 'Evaluate should match manual computation');

%% 3) Reach and soundness test (node-only perturbation)
GS_out = L.reach(GS_in, E, adj_list, 'approx-star');
assert(isa(GS_out, 'GraphStar'), 'Output should be GraphStar');
assert(GS_out.numNodes == numNodes, 'Output should have same number of nodes');
assert(GS_out.numFeatures == 8, 'Output should have 8 features');

% Soundness: center should match evaluate
center_in = GS_in.V(:, :, 1);
center_out = GS_out.V(:, :, 1);
expected_center = L.evaluate(center_in, E, adj_list);
assert(max(abs(center_out - expected_center), [], 'all') < 1e-10, 'Center should match evaluate');

% Containment: center output should be within bounds
[lb_out, ub_out] = GS_out.getRanges();
Y_center = L.evaluate(GS_in.V(:,:,1), E, adj_list);
tol = 1e-6;
assert(all(Y_center(:) >= lb_out(:) - tol), 'Center output should be >= lower bound');
assert(all(Y_center(:) <= ub_out(:) + tol), 'Center output should be <= upper bound');

%% 4) Edge perturbation reach and soundness test
E_lb = E - 0.05;
E_ub = E + 0.05;
E_star = ImageStar(E_lb, E_ub);

GS_out_edge = L.reach(GS_in, E_star, adj_list, 'approx-star');
assert(isa(GS_out_edge, 'GraphStar'), 'Output should be GraphStar with edge perturbation');
assert(GS_out_edge.numNodes == numNodes, 'Output should have same number of nodes');

% Soundness: center should match evaluate
center_in_edge = GS_in.V(:, :, 1);  % Recompute center_in for this section
center_out_edge = GS_out_edge.V(:, :, 1);
expected_center_edge = L.evaluate(center_in_edge, E, adj_list);
assert(max(abs(center_out_edge - expected_center_edge), [], 'all') < 1e-10, ...
    'Center should match evaluate for edge perturbation');

% Containment: center output should be within bounds
[lb_out_edge, ub_out_edge] = GS_out_edge.getRanges();
Y_center_edge = L.evaluate(GS_in.V(:,:,1), E, adj_list);
tol_edge = 1e-6;
assert(all(Y_center_edge(:) >= lb_out_edge(:) - tol_edge), 'Center output should be >= lower bound (edge perturbation)');
assert(all(Y_center_edge(:) <= ub_out_edge(:) + tol_edge), 'Center output should be <= upper bound (edge perturbation)');

%% 5) Precision change test
L_prec = GINELayer('prec_test', W_node, b_node, W_edge, b_edge);
L_prec.changeParamsPrecision('single');
assert(isa(L_prec.Weights, 'single'), 'Node weights should be single precision');
assert(isa(L_prec.EdgeWeights, 'single'), 'Edge weights should be single precision');

L_prec.changeParamsPrecision('double');
assert(isa(L_prec.Weights, 'double'), 'Node weights should be double precision');

disp('All GINELayer tests passed!');
