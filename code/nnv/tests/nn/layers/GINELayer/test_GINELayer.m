% test GINELayer class

%% 1) Constructor test - 5 arguments (name, W_node, b_node, W_edge, b_edge)
W_node = rand(4, 8);     % 4 input -> 8 output features
b_node = rand(8, 1);
W_edge = rand(3, 4);     % 3 edge features -> 4 (node input dim)
b_edge = rand(4, 1);
L = GINELayer('test_gine', W_node, b_node, W_edge, b_edge);

assert(L.InputSize == 4, 'InputSize should be 4');
assert(L.OutputSize == 8, 'OutputSize should be 8');
assert(L.EdgeInputSize == 3, 'EdgeInputSize should be 3');
assert(strcmp(L.Name, 'test_gine'), 'Name should match');
assert(isequal(L.Weights, W_node), 'Node weights should match');
assert(isequal(L.EdgeWeights, W_edge), 'Edge weights should match');
assert(L.Epsilon == 0, 'Default epsilon should be 0');

%% 2) Constructor test - 4 arguments (W_node, b_node, W_edge, b_edge)
L2 = GINELayer(W_node, b_node, W_edge, b_edge);

assert(L2.InputSize == 4, 'InputSize should be 4');
assert(L2.OutputSize == 8, 'OutputSize should be 8');
assert(strcmp(L2.Name, 'gine_layer'), 'Default name should be gine_layer');

%% 3) Constructor test - 6 arguments (with epsilon)
L3 = GINELayer('eps_test', W_node, b_node, W_edge, b_edge, 0.1);

assert(L3.Epsilon == 0.1, 'Epsilon should be 0.1');

%% 4) Constructor test - 0 arguments (empty)
L0 = GINELayer();

assert(L0.InputSize == 0, 'Empty layer InputSize should be 0');
assert(L0.OutputSize == 0, 'Empty layer OutputSize should be 0');
assert(L0.EdgeInputSize == 0, 'Empty layer EdgeInputSize should be 0');
assert(isempty(L0.Weights), 'Empty layer Weights should be empty');

%% 5) Evaluate test with simple graph
numNodes = 4;
numEdges = 5;
X = rand(numNodes, 4);          % node features
E = rand(numEdges, 3);          % edge features
adj_list = [1 2; 1 3; 2 3; 3 4; 4 1];  % edges

Y = L.evaluate(X, E, adj_list);

assert(size(Y, 1) == numNodes, 'Output should have same number of nodes');
assert(size(Y, 2) == 8, 'Output should have 8 features');

%% 6) Verify evaluate computation manually
% Compute expected output step by step
E_trans = E * W_edge + b_edge';
src_nodes = adj_list(:, 1);
X_src = X(src_nodes, :);
edge_msg = max(0, X_src + E_trans);  % ReLU

agg = zeros(numNodes, 4);
dst_nodes = adj_list(:, 2);
for e = 1:numEdges
    agg(dst_nodes(e), :) = agg(dst_nodes(e), :) + edge_msg(e, :);
end

combined = (1 + L.Epsilon) * X + agg;
expected_Y = combined * W_node + b_node';

assert(max(abs(Y - expected_Y), [], 'all') < 1e-10, 'Evaluate should match manual computation');

%% 7) GraphStar reachability test - node-only mode
NF = rand(numNodes, 4);
LB = -0.1 * ones(numNodes, 4);
UB = 0.1 * ones(numNodes, 4);

GS_in = GraphStar(NF, LB, UB);

GS_out = L.reach(GS_in, E, adj_list, 'approx-star');

assert(isa(GS_out, 'GraphStar'), 'Output should be GraphStar');
assert(GS_out.numNodes == numNodes, 'Output should have same number of nodes');
assert(GS_out.numFeatures == 8, 'Output should have 8 features');

%% 8) Verify center matches evaluate
center_in = GS_in.V(:, :, 1);  % center of input GraphStar
center_out = GS_out.V(:, :, 1);  % center of output GraphStar
expected_center = L.evaluate(center_in, E, adj_list);

assert(max(abs(center_out - expected_center), [], 'all') < 1e-10, 'Center should match evaluate');

%% 9) Test without explicit method (should default to approx-star)
GS_out2 = L.reach(GS_in, E, adj_list);

assert(isa(GS_out2, 'GraphStar'), 'Output should be GraphStar with default method');

%% 10) Test with different epsilon value
L_eps = GINELayer('eps_layer', W_node, b_node, W_edge, b_edge, 0.5);
Y_eps = L_eps.evaluate(X, E, adj_list);

% With epsilon=0.5, self-loop contribution is 1.5*X instead of 1.0*X
combined_eps = 1.5 * X + agg;
expected_Y_eps = combined_eps * W_node + b_node';

assert(max(abs(Y_eps - expected_Y_eps), [], 'all') < 1e-10, 'Epsilon scaling should work');

%% 11) Test precision change
L_single = GINELayer('single_test', W_node, b_node, W_edge, b_edge);
L_single.changeParamsPrecision('single');

assert(isa(L_single.Weights, 'single'), 'Node weights should be single precision');
assert(isa(L_single.EdgeWeights, 'single'), 'Edge weights should be single precision');
assert(isa(L_single.Bias, 'single'), 'Node bias should be single precision');
assert(isa(L_single.EdgeBias, 'single'), 'Edge bias should be single precision');

L_single.changeParamsPrecision('double');
assert(isa(L_single.Weights, 'double'), 'Node weights should be double precision');

%% 12) Test output bounds contain samples
% Generate random samples from input and verify they're in output bounds
num_samples = 10;
for s = 1:num_samples
    % Generate random sample in input bounds
    alpha = rand(GS_in.numPred, 1) .* (GS_in.pred_ub - GS_in.pred_lb) + GS_in.pred_lb;
    X_sample = GS_in.V(:, :, 1);
    for k = 1:GS_in.numPred
        X_sample = X_sample + alpha(k) * GS_in.V(:, :, k+1);
    end

    % Evaluate at sample
    Y_sample = L.evaluate(X_sample, E, adj_list);

    % Get output bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Check sample is within bounds (with small tolerance for numerical error)
    tol = 1e-6;
    assert(all(Y_sample(:) >= lb_out(:) - tol), 'Sample should be above lower bound');
    assert(all(Y_sample(:) <= ub_out(:) + tol), 'Sample should be below upper bound');
end

%% ========== Phase 2: Edge Perturbation Tests ==========

%% 13) Edge perturbation mode with ImageStar
% Create edge perturbation as ImageStar
E_center = E;
E_lb = E - 0.05;
E_ub = E + 0.05;
E_star = ImageStar(E_lb, E_ub);

GS_out_edge = L.reach(GS_in, E_star, adj_list, 'approx-star');

assert(isa(GS_out_edge, 'GraphStar'), 'Output should be GraphStar with edge perturbation');
assert(GS_out_edge.numNodes == numNodes, 'Output should have same number of nodes');
assert(GS_out_edge.numFeatures == 8, 'Output should have 8 features');

%% 14) Verify center matches evaluate for edge perturbation
center_in_node = GS_in.V(:, :, 1);  % center of input node GraphStar
center_out_edge = GS_out_edge.V(:, :, 1);  % center of output GraphStar
expected_center_edge = L.evaluate(center_in_node, E_center, adj_list);

assert(max(abs(center_out_edge - expected_center_edge), [], 'all') < 1e-10, ...
    'Center should match evaluate for edge perturbation');

%% 15) Test that edge perturbation produces larger bounds than node-only
% With edge perturbation, bounds should generally be larger
[lb_node_only, ub_node_only] = GS_out.getRanges();
[lb_edge_pert, ub_edge_pert] = GS_out_edge.getRanges();

% The edge perturbation bounds should contain the node-only bounds
% (allowing small tolerance for numerical differences)
tol = 1e-6;
assert(all(lb_edge_pert(:) <= lb_node_only(:) + tol), ...
    'Edge perturbation lower bounds should be at most node-only bounds');
assert(all(ub_edge_pert(:) >= ub_node_only(:) - tol), ...
    'Edge perturbation upper bounds should be at least node-only bounds');

%% 16) Test edge perturbation with samples
num_samples_edge = 10;
for s = 1:num_samples_edge
    % Generate random sample for node features
    alpha_node = rand(GS_in.numPred, 1) .* (GS_in.pred_ub - GS_in.pred_lb) + GS_in.pred_lb;
    X_sample = GS_in.V(:, :, 1);
    for k = 1:GS_in.numPred
        X_sample = X_sample + alpha_node(k) * GS_in.V(:, :, k+1);
    end

    % Generate random sample for edge features (within perturbation bounds)
    E_sample = E_lb + rand(size(E)) .* (E_ub - E_lb);

    % Evaluate at sample
    Y_sample = L.evaluate(X_sample, E_sample, adj_list);

    % Check sample is within bounds (with tolerance for over-approximation)
    tol = 1e-5;
    assert(all(Y_sample(:) >= lb_edge_pert(:) - tol), ...
        'Sample should be above edge perturbation lower bound');
    assert(all(Y_sample(:) <= ub_edge_pert(:) + tol), ...
        'Sample should be below edge perturbation upper bound');
end

%% 17) Test edge perturbation preserves output type
assert(isa(GS_out_edge, 'GraphStar'), 'Edge perturbation output should be GraphStar');
assert(GS_out_edge.numPred >= GS_in.numPred, ...
    'Edge perturbation should have at least as many predicates as input');

disp('All GINELayer tests passed (including Phase 2 edge perturbation)!');
