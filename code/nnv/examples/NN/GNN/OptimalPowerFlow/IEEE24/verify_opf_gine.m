% verify_opf_gine.m - GINE Optimal Power Flow Verification Example (IEEE 24-bus)
%
% This example demonstrates GNN verification for Optimal Power Flow prediction
% using a 3-layer GINE model on the IEEE 24-bus system.
%
% The script performs:
%   1. Reachability analysis to compute output bounds
%   2. Voltage magnitude specification verification (0.95-1.05 p.u.)
%
% Author: Anne Tumlin
% Date: 01/13/2026

%% Setup
clear; clc;

% Load the trained GINE model
modelPath = fullfile(fileparts(mfilename('fullpath')), 'models', 'gine_opf_ieee24.mat');
model = load(modelPath);

fprintf('=== GINE Optimal Power Flow Verification (IEEE 24-bus) ===\n');
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
edge_weights = double(model.a);

X = double(model.X_test_g{1});
numNodes = size(X, 1);
numEdges = size(adj_list, 1);

fprintf('\nGraph structure: %d nodes, %d edges\n', numNodes, numEdges);

%% Create GNN wrapper
gnn = GNN({L1, L2, L3}, [], adj_list, E, edge_weights);

%% Perturbation configuration (matching GNNV prototype)
% Which features to perturb:
%   [] = all features (default)
%   [1, 2] = only power injections (Pg-Pd, Qg-Qd) - matches GNNV config
perturb_features = [1, 2];  % Match GNNV: perturb only power injections
epsilon = 0.01;
range_per_col = max(X) - min(X);

% Define input perturbation (selective features)
eps_matrix = zeros(numNodes, size(X, 2));
if isempty(perturb_features)
    scaled_eps = range_per_col .* epsilon;
    eps_matrix = repmat(scaled_eps, numNodes, 1);
else
    for f = perturb_features
        if f <= size(X, 2)
            eps_matrix(:, f) = range_per_col(f) * epsilon;
        end
    end
end

%% Create GraphStar and compute reachability
% GraphStar(NF, LB, UB) expects perturbation bounds relative to NF
% So we pass -eps_matrix and +eps_matrix (not absolute bounds)
GS_in = GraphStar(X, -eps_matrix, eps_matrix);

fprintf('\n=== Computing Reachability ===\n');
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
    alpha = rand(GS_in.numPred, 1) .* (GS_in.pred_ub - GS_in.pred_lb) + GS_in.pred_lb;
    X_sample = GS_in.V(:, :, 1);
    for k = 1:GS_in.numPred
        X_sample = X_sample + alpha(k) * GS_in.V(:, :, k+1);
    end
    Y_sample = gnn.evaluate(X_sample);
    tol = 1e-6;
    if all(Y_sample(:) >= lb_out(:) - tol) && all(Y_sample(:) <= ub_out(:) + tol)
        samples_in_bounds = samples_in_bounds + 1;
    end
end
fprintf('Samples within bounds: %d/%d\n', samples_in_bounds, num_samples);

%% Verify voltage magnitude specification
v_min = 0.95;  % Per-unit lower bound
v_max = 1.05;  % Per-unit upper bound

% Add PowerFlow folder to path for verify_voltage_spec
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'PowerFlow'));
results = verify_voltage_spec(GS_out, model, v_min, v_max);

fprintf('\n=== Voltage Specification Verification ===\n');
fprintf('Specification: %.2f <= V <= %.2f p.u.\n', v_min, v_max);
fprintf('  Verified safe: %d nodes\n', sum(results == 1));
fprintf('  Violated: %d nodes\n', sum(results == 0));
unknown_boundary = sum(results == 2);
unknown_timeout = sum(results == 3);
fprintf('  Unknown: %d nodes\n', unknown_boundary + unknown_timeout);
if unknown_boundary > 0 || unknown_timeout > 0
    fprintf('    - Bounds cross spec boundary: %d\n', unknown_boundary);
    fprintf('    - Timeout/inconclusive: %d\n', unknown_timeout);
end
fprintf('  N/A (non-voltage bus): %d nodes\n', sum(results == -1));

fprintf('\n=== Complete ===\n');
