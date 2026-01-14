% verify_pf_gine_edge_perturb.m - GINE Power Flow Verification with Edge Perturbation (IEEE 39-bus)
%
% This example demonstrates GNN verification for Power Flow prediction
% using a 3-layer GINE model on the IEEE 39-bus (New England) system with
% BOTH node and edge feature perturbations.
%
% The script performs:
%   1. Reachability analysis with node AND edge perturbations
%   2. Voltage magnitude specification verification (0.95-1.05 p.u.)
%
% This matches the GNNV prototype edge perturbation configuration.
%
% Author: Anne Tumlin
% Date: 01/14/2026

%% Setup
clear; clc;

% Load the trained GINE model
modelPath = fullfile(fileparts(mfilename('fullpath')), 'models', 'gine_edgelist_pf_ieee39_run3_seed130.mat');
model = load(modelPath);

fprintf('=== GINE Power Flow Verification with Edge Perturbation (IEEE 39-bus) ===\n');
fprintf('Model: %s\n\n', modelPath);

%% Extract model weights and create GINE layers
params = model.best_params;

W_node1 = double(gather(params.mult1.Weights));
W_node2 = double(gather(params.mult2.Weights));
W_node3 = double(gather(params.mult3.Weights));

W_edge1 = double(gather(params.edge1.Weights));
W_edge2 = double(gather(params.edge2.Weights));
W_edge3 = double(gather(params.edge3.Weights));

b_node1 = zeros(size(W_node1, 2), 1);
b_node2 = zeros(size(W_node2, 2), 1);
b_node3 = zeros(size(W_node3, 2), 1);
b_edge1 = zeros(size(W_edge1, 2), 1);
b_edge2 = zeros(size(W_edge2, 2), 1);
b_edge3 = zeros(size(W_edge3, 2), 1);

fprintf('Layer dimensions:\n');
fprintf('  Layer 1: Node %d->%d, Edge %d->%d\n', size(W_node1, 1), size(W_node1, 2), size(W_edge1, 1), size(W_edge1, 2));
fprintf('  Layer 2: Node %d->%d, Edge %d->%d\n', size(W_node2, 1), size(W_node2, 2), size(W_edge2, 1), size(W_edge2, 2));
fprintf('  Layer 3: Node %d->%d, Edge %d->%d\n', size(W_node3, 1), size(W_node3, 2), size(W_edge3, 1), size(W_edge3, 2));

L1 = GINELayer('gine1', W_node1, b_node1, W_edge1, b_edge1);
L2 = GINELayer('gine2', W_node2, b_node2, W_edge2, b_edge2);
L3 = GINELayer('gine3', W_node3, b_node3, W_edge3, b_edge3);

%% Extract graph structure
src = double(model.src);
dst = double(model.dst);
adj_list = [src, dst];
E = double(model.E_edge);
edge_weights = double(model.a);  % Normalized adjacency weights

X = double(model.X_test_g{1});
numNodes = size(X, 1);
numEdges = size(adj_list, 1);

fprintf('\nGraph structure: %d nodes, %d edges\n', numNodes, numEdges);
fprintf('Edge features: %d edges x %d features\n', size(E, 1), size(E, 2));

%% Perturbation configuration
% Node feature perturbation settings
epsilon_node = 0.001;  % 0.1% perturbation
perturb_node_features = [1, 2];  % Power injections (Pg-Pd, Qg-Qd)

% Edge feature perturbation settings
epsilon_edge = 0.001;  % 0.1% perturbation
perturb_edge_features = [1];  % First edge feature (impedance)

fprintf('\n=== Perturbation Configuration ===\n');
fprintf('Node epsilon: %.3f, features: %s\n', epsilon_node, mat2str(perturb_node_features));
fprintf('Edge epsilon: %.3f, features: %s\n', epsilon_edge, mat2str(perturb_edge_features));

%% Create Node Perturbation (GraphStar)
range_per_col = max(X) - min(X);
eps_matrix_node = zeros(numNodes, size(X, 2));
for f = perturb_node_features
    if f <= size(X, 2)
        eps_matrix_node(:, f) = range_per_col(f) * epsilon_node;
    end
end
GS_in = GraphStar(X, -eps_matrix_node, eps_matrix_node);

%% Create Edge Perturbation (GraphStar - consistent API with nodes)
range_per_edge_col = max(E) - min(E);
eps_matrix_edge = zeros(numEdges, size(E, 2));
for f = perturb_edge_features
    if f <= size(E, 2)
        eps_matrix_edge(:, f) = range_per_edge_col(f) * epsilon_edge;
    end
end
E_star = GraphStar(E, -eps_matrix_edge, eps_matrix_edge);

fprintf('\nNode perturbation: %d nodes x %d features\n', size(eps_matrix_node, 1), size(eps_matrix_node, 2));
fprintf('Edge perturbation: %d edges x %d features\n', size(eps_matrix_edge, 1), size(eps_matrix_edge, 2));

%% Create GNN with Edge Star (enables edge perturbation mode)
gnn = GNN({L1, L2, L3}, [], adj_list, E_star, edge_weights);

%% Compute Reachability
fprintf('\n=== Computing Reachability (Node + Edge Perturbation) ===\n');
reachOpts = struct('reachMethod', 'approx-star', 'dis_opt', 'display');
t_start = tic;
GS_out = gnn.reach(GS_in, reachOpts);
total_time = toc(t_start);

fprintf('\nCompleted in %.4f seconds\n', total_time);

%% Analyze output bounds
[lb_out, ub_out] = GS_out.getRanges();
bound_widths = ub_out - lb_out;

fprintf('\nOutput bounds - Mean: %.6f, Max: %.6f\n', mean(bound_widths(:)), max(bound_widths(:)));

%% Sample validation
num_samples = 10;
samples_in_bounds = 0;
for s = 1:num_samples
    % Sample node features
    alpha_node = rand(GS_in.numPred, 1) .* (GS_in.pred_ub - GS_in.pred_lb) + GS_in.pred_lb;
    X_sample = GS_in.V(:, :, 1);
    for k = 1:GS_in.numPred
        X_sample = X_sample + alpha_node(k) * GS_in.V(:, :, k+1);
    end

    % Sample edge features
    alpha_edge = rand(E_star.numPred, 1) .* (E_star.pred_ub - E_star.pred_lb) + E_star.pred_lb;
    E_sample = E_star.V(:, :, 1);
    for k = 1:E_star.numPred
        E_sample = E_sample + alpha_edge(k) * E_star.V(:, :, k+1);
    end

    % Evaluate with sampled inputs
    Y_sample = gnn.evaluate(X_sample, E_sample);

    tol = 1e-5;
    if all(Y_sample(:) >= lb_out(:) - tol) && all(Y_sample(:) <= ub_out(:) + tol)
        samples_in_bounds = samples_in_bounds + 1;
    end
end
fprintf('Samples within bounds: %d/%d\n', samples_in_bounds, num_samples);

%% Verify voltage magnitude specification
v_min = 0.95;  % Per-unit lower bound
v_max = 1.05;  % Per-unit upper bound

% Add parent folder to path for verify_voltage_spec
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
results = verify_voltage_spec(GS_out, model, v_min, v_max);

fprintf('\n=== Voltage Specification Verification ===\n');
fprintf('Specification: %.2f <= V <= %.2f p.u.\n', v_min, v_max);
fprintf('  Verified safe: %d nodes\n', sum(results == 1));
fprintf('  Violated: %d nodes\n', sum(results == 0));
fprintf('  Unknown: %d nodes\n', sum(results == 2));
fprintf('  N/A (non-voltage bus): %d nodes\n', sum(results == -1));

fprintf('\n=== Complete (Edge Perturbation Mode) ===\n');
