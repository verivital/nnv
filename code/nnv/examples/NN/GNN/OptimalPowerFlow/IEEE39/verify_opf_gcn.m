% verify_opf_gcn.m - GCN Optimal Power Flow Verification Example (IEEE 39-bus)
%
% This example demonstrates GNN verification for Optimal Power Flow prediction
% using a 3-layer GCN model on the IEEE 39-bus (New England) system.
%
% The script performs:
%   1. Reachability analysis to compute output bounds
%   2. Voltage magnitude specification verification (0.95-1.05 p.u.)
%
% Author: Anne Tumlin
% Date: 01/13/2026

%% Setup
clear; clc;

% Load the trained GCN model
modelPath = fullfile(fileparts(mfilename('fullpath')), 'models', 'gcn_opf_ieee39_run1_seed128.mat');
model = load(modelPath);

fprintf('=== GCN Optimal Power Flow Verification (IEEE 39-bus) ===\n');
fprintf('Model: %s\n\n', modelPath);

%% Extract model weights and create GCN layers
params = model.best_params;

W1 = double(gather(params.mult1.Weights));
W2 = double(gather(params.mult2.Weights));
W3 = double(gather(params.mult3.Weights));

b1 = zeros(size(W1, 2), 1);
b2 = zeros(size(W2, 2), 1);
b3 = zeros(size(W3, 2), 1);

fprintf('Layer dimensions:\n');
fprintf('  Layer 1: %d -> %d\n', size(W1, 1), size(W1, 2));
fprintf('  Layer 2: %d -> %d\n', size(W2, 1), size(W2, 2));
fprintf('  Layer 3: %d -> %d\n', size(W3, 1), size(W3, 2));

L1 = GCNLayer('gcn1', W1, b1);
R1 = ReluLayer();  % ReLU after layer 1
L2 = GCNLayer('gcn2', W2, b2);
R2 = ReluLayer();  % ReLU after layer 2
L3 = GCNLayer('gcn3', W3, b3);
R3 = ReluLayer();  % ReLU after layer 3 (matches GNNV)

%% Extract graph structure
A_norm = double(model.ANorm_g);
numNodes = size(A_norm, 1);

fprintf('\nGraph structure: %d nodes\n', numNodes);

%% Create GNN wrapper (with ReLU activations)
gnn = GNN({L1, R1, L2, R2, L3, R3}, A_norm);

%% Prepare test data
X = double(model.X_test_g{1});

fprintf('Input: %d nodes x %d features\n', size(X, 1), size(X, 2));

%% Define input perturbation (1% of feature range)
epsilon = 0.01;
range_per_col = max(X) - min(X);
scaled_eps = range_per_col .* epsilon;
eps_matrix = repmat(scaled_eps, numNodes, 1);

lb = X - eps_matrix;
ub = X + eps_matrix;

%% Create GraphStar and compute reachability
GS_in = GraphStar(X, lb, ub);

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
fprintf('  Unknown: %d nodes\n', sum(results == 2));
fprintf('  N/A (non-voltage bus): %d nodes\n', sum(results == -1));

fprintf('\n=== Complete ===\n');
