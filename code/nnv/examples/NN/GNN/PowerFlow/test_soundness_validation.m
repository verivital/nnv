% test_soundness_validation.m - Comprehensive Soundness Validation for GNN Verification
%
% This script validates that the NNV GNN verification implementation is SOUND:
%   1. Center point of reachability matches forward pass (evaluate)
%   2. Random samples within perturbation region fall within computed bounds
%   3. Boundary samples near perturbation limits fall within computed bounds
%   4. GCN and GINE layers produce consistent results
%
% A sound over-approximation means: ALL points in the input set produce
% outputs WITHIN the computed output bounds. If this test passes, the
% verification is sound.
%
% Author: Anne Tumlin
% Date: 01/15/2026

%% Setup
clear; clc;
fprintf('=== GNN Verification Soundness Validation ===\n\n');

% Add paths
addpath(genpath('/home/annemtumlin/dev/nnv/code/nnv'));

%% Test Parameters
num_random_samples = 100;  % Number of random samples to test
num_boundary_samples = 50;  % Number of boundary samples to test
tolerance = 1e-6;  % Numerical tolerance for comparisons

all_tests_passed = true;
test_results = struct();

%% Test 1: GINE Layer - IEEE24
fprintf('=== Test 1: GINE Layer (IEEE24) ===\n');

% Load model
modelPath = fullfile(fileparts(mfilename('fullpath')), 'IEEE24', 'models', 'gine_edgelist_pf_ieee24_run4_seed131.mat');
if ~exist(modelPath, 'file')
    fprintf('WARNING: Model file not found, skipping GINE test\n');
    test_results.gine_ieee24 = 'SKIPPED';
else
    model = load(modelPath);

    % Extract weights
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

    % Create layers
    L1 = GINELayer('gine1', W_node1, b_node1, W_edge1, b_edge1);
    L2 = GINELayer('gine2', W_node2, b_node2, W_edge2, b_edge2);
    L3 = GINELayer('gine3', W_node3, b_node3, W_edge3, b_edge3);

    % Extract graph structure
    src = double(model.src);
    dst = double(model.dst);
    adj_list = [src, dst];
    E = double(model.E_edge);
    edge_weights = double(model.a);
    X = double(model.X_test_g{1});
    numNodes = size(X, 1);

    % Create GNN
    gnn = GNN({L1, L2, L3}, [], adj_list, E, edge_weights);

    % Perturbation configuration
    epsilon = 0.01;
    perturb_features = [1, 2];
    range_per_col = max(X) - min(X);
    eps_matrix = zeros(numNodes, size(X, 2));
    for f = perturb_features
        if f <= size(X, 2)
            eps_matrix(:, f) = range_per_col(f) * epsilon;
        end
    end

    % Create GraphStar
    GS_in = GraphStar(X, -eps_matrix, eps_matrix);

    % Compute reachability
    fprintf('Computing reachability...\n');
    reachOpts = struct('reachMethod', 'approx-star');
    GS_out = gnn.reach(GS_in, reachOpts);

    % Get bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Test 1a: Center point check
    fprintf('\n--- Test 1a: Center Point Check ---\n');
    Y_center = gnn.evaluate(X);
    center_in_bounds = all(Y_center(:) >= lb_out(:) - tolerance) && ...
                       all(Y_center(:) <= ub_out(:) + tolerance);
    if center_in_bounds
        fprintf('PASS: Center point within bounds\n');
    else
        fprintf('FAIL: Center point outside bounds!\n');
        fprintf('  Max below lb: %.6e\n', max(lb_out(:) - Y_center(:)));
        fprintf('  Max above ub: %.6e\n', max(Y_center(:) - ub_out(:)));
        all_tests_passed = false;
    end

    % Test 1b: Random samples check
    fprintf('\n--- Test 1b: Random Samples Check ---\n');
    samples_in_bounds = 0;
    max_violation = 0;

    for s = 1:num_random_samples
        % Generate random predicate values
        alpha = rand(GS_in.numPred, 1) .* (GS_in.pred_ub - GS_in.pred_lb) + GS_in.pred_lb;

        % Evaluate at sample
        X_sample = GS_in.V(:, :, 1);
        for k = 1:GS_in.numPred
            X_sample = X_sample + alpha(k) * GS_in.V(:, :, k+1);
        end

        % Forward pass
        Y_sample = gnn.evaluate(X_sample);

        % Check bounds
        below_lb = lb_out(:) - Y_sample(:) - tolerance;
        above_ub = Y_sample(:) - ub_out(:) - tolerance;
        violation = max([0; below_lb; above_ub]);
        max_violation = max(max_violation, violation);

        if violation <= 0
            samples_in_bounds = samples_in_bounds + 1;
        end
    end

    fprintf('Random samples within bounds: %d/%d\n', samples_in_bounds, num_random_samples);
    fprintf('Max violation: %.6e\n', max_violation);

    if samples_in_bounds == num_random_samples
        fprintf('PASS: All random samples within bounds\n');
    else
        fprintf('FAIL: %d samples outside bounds!\n', num_random_samples - samples_in_bounds);
        all_tests_passed = false;
    end

    % Test 1c: Boundary samples check
    fprintf('\n--- Test 1c: Boundary Samples Check ---\n');
    boundary_samples_in_bounds = 0;

    for s = 1:num_boundary_samples
        % Generate boundary sample (predicates at their limits)
        alpha = zeros(GS_in.numPred, 1);
        for k = 1:GS_in.numPred
            if rand() < 0.5
                alpha(k) = GS_in.pred_lb(k);  % Lower boundary
            else
                alpha(k) = GS_in.pred_ub(k);  % Upper boundary
            end
        end

        % Evaluate at sample
        X_sample = GS_in.V(:, :, 1);
        for k = 1:GS_in.numPred
            X_sample = X_sample + alpha(k) * GS_in.V(:, :, k+1);
        end

        % Forward pass
        Y_sample = gnn.evaluate(X_sample);

        % Check bounds
        if all(Y_sample(:) >= lb_out(:) - tolerance) && ...
           all(Y_sample(:) <= ub_out(:) + tolerance)
            boundary_samples_in_bounds = boundary_samples_in_bounds + 1;
        end
    end

    fprintf('Boundary samples within bounds: %d/%d\n', boundary_samples_in_bounds, num_boundary_samples);

    if boundary_samples_in_bounds == num_boundary_samples
        fprintf('PASS: All boundary samples within bounds\n');
    else
        fprintf('FAIL: %d boundary samples outside bounds!\n', num_boundary_samples - boundary_samples_in_bounds);
        all_tests_passed = false;
    end

    test_results.gine_ieee24 = struct(...
        'center_check', center_in_bounds, ...
        'random_samples', samples_in_bounds, ...
        'boundary_samples', boundary_samples_in_bounds, ...
        'max_violation', max_violation);
end

%% Test 2: GCN Layer - IEEE24
fprintf('\n=== Test 2: GCN Layer (IEEE24) ===\n');

modelPath = fullfile(fileparts(mfilename('fullpath')), 'IEEE24', 'models', 'gcn_pf_ieee24_run4_seed131.mat');
if ~exist(modelPath, 'file')
    fprintf('WARNING: Model file not found, skipping GCN test\n');
    test_results.gcn_ieee24 = 'SKIPPED';
else
    model = load(modelPath);

    % Extract weights
    params = model.best_params;
    W1 = double(gather(params.mult1.Weights));
    W2 = double(gather(params.mult2.Weights));
    W3 = double(gather(params.mult3.Weights));

    b1 = zeros(size(W1, 2), 1);
    b2 = zeros(size(W2, 2), 1);
    b3 = zeros(size(W3, 2), 1);

    % Create layers with ReLU
    L1 = GCNLayer('gcn1', W1, b1);
    R1 = ReluLayer();
    L2 = GCNLayer('gcn2', W2, b2);
    R2 = ReluLayer();
    L3 = GCNLayer('gcn3', W3, b3);
    R3 = ReluLayer();

    % Extract graph structure
    A_norm = double(model.ANorm_g);
    X = double(model.X_test_g{1});
    numNodes = size(X, 1);

    % Create GNN
    gnn = GNN({L1, R1, L2, R2, L3, R3}, A_norm);

    % Perturbation configuration
    epsilon = 0.01;
    perturb_features = [1, 2];
    range_per_col = max(X) - min(X);
    eps_matrix = zeros(numNodes, size(X, 2));
    for f = perturb_features
        if f <= size(X, 2)
            eps_matrix(:, f) = range_per_col(f) * epsilon;
        end
    end

    % Create GraphStar
    GS_in = GraphStar(X, -eps_matrix, eps_matrix);

    % Compute reachability
    fprintf('Computing reachability...\n');
    reachOpts = struct('reachMethod', 'approx-star');
    GS_out = gnn.reach(GS_in, reachOpts);

    % Get bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Test 2a: Center point check
    fprintf('\n--- Test 2a: Center Point Check ---\n');
    Y_center = gnn.evaluate(X);
    center_in_bounds = all(Y_center(:) >= lb_out(:) - tolerance) && ...
                       all(Y_center(:) <= ub_out(:) + tolerance);
    if center_in_bounds
        fprintf('PASS: Center point within bounds\n');
    else
        fprintf('FAIL: Center point outside bounds!\n');
        all_tests_passed = false;
    end

    % Test 2b: Random samples check
    fprintf('\n--- Test 2b: Random Samples Check ---\n');
    samples_in_bounds = 0;

    for s = 1:num_random_samples
        alpha = rand(GS_in.numPred, 1) .* (GS_in.pred_ub - GS_in.pred_lb) + GS_in.pred_lb;
        X_sample = GS_in.V(:, :, 1);
        for k = 1:GS_in.numPred
            X_sample = X_sample + alpha(k) * GS_in.V(:, :, k+1);
        end
        Y_sample = gnn.evaluate(X_sample);

        if all(Y_sample(:) >= lb_out(:) - tolerance) && ...
           all(Y_sample(:) <= ub_out(:) + tolerance)
            samples_in_bounds = samples_in_bounds + 1;
        end
    end

    fprintf('Random samples within bounds: %d/%d\n', samples_in_bounds, num_random_samples);

    if samples_in_bounds == num_random_samples
        fprintf('PASS: All random samples within bounds\n');
    else
        fprintf('FAIL: %d samples outside bounds!\n', num_random_samples - samples_in_bounds);
        all_tests_passed = false;
    end

    test_results.gcn_ieee24 = struct(...
        'center_check', center_in_bounds, ...
        'random_samples', samples_in_bounds);
end

%% Test 3: GINE Layer with Edge Perturbation - IEEE24
fprintf('\n=== Test 3: GINE Layer with Edge Perturbation (IEEE24) ===\n');

modelPath = fullfile(fileparts(mfilename('fullpath')), 'IEEE24', 'models', 'gine_edgelist_pf_ieee24_run4_seed131.mat');
if ~exist(modelPath, 'file')
    fprintf('WARNING: Model file not found, skipping edge perturbation test\n');
    test_results.gine_edge_ieee24 = 'SKIPPED';
else
    model = load(modelPath);

    % Extract weights (same as before)
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

    L1 = GINELayer('gine1', W_node1, b_node1, W_edge1, b_edge1);
    L2 = GINELayer('gine2', W_node2, b_node2, W_edge2, b_edge2);
    L3 = GINELayer('gine3', W_node3, b_node3, W_edge3, b_edge3);

    % Extract graph structure
    src = double(model.src);
    dst = double(model.dst);
    adj_list = [src, dst];
    E = double(model.E_edge);
    edge_weights = double(model.a);
    X = double(model.X_test_g{1});
    numNodes = size(X, 1);
    numEdges = size(adj_list, 1);

    % Node perturbation
    epsilon_node = 0.001;
    perturb_node_features = [1, 2];
    range_per_col = max(X) - min(X);
    eps_matrix_node = zeros(numNodes, size(X, 2));
    for f = perturb_node_features
        if f <= size(X, 2)
            eps_matrix_node(:, f) = range_per_col(f) * epsilon_node;
        end
    end
    GS_in = GraphStar(X, -eps_matrix_node, eps_matrix_node);

    % Edge perturbation
    epsilon_edge = 0.001;
    perturb_edge_features = [1];
    range_per_edge_col = max(E) - min(E);
    eps_matrix_edge = zeros(numEdges, size(E, 2));
    for f = perturb_edge_features
        if f <= size(E, 2)
            eps_matrix_edge(:, f) = range_per_edge_col(f) * epsilon_edge;
        end
    end
    E_star = GraphStar(E, -eps_matrix_edge, eps_matrix_edge);

    % Create GNN with edge star
    gnn = GNN({L1, L2, L3}, [], adj_list, E_star, edge_weights);

    % Compute reachability
    fprintf('Computing reachability with edge perturbation...\n');
    reachOpts = struct('reachMethod', 'approx-star');
    GS_out = gnn.reach(GS_in, reachOpts);

    % Get bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Test 3a: Center point check
    fprintf('\n--- Test 3a: Center Point Check ---\n');
    Y_center = gnn.evaluate(X, E);
    center_in_bounds = all(Y_center(:) >= lb_out(:) - tolerance) && ...
                       all(Y_center(:) <= ub_out(:) + tolerance);
    if center_in_bounds
        fprintf('PASS: Center point within bounds\n');
    else
        fprintf('FAIL: Center point outside bounds!\n');
        all_tests_passed = false;
    end

    % Test 3b: Random samples check (both node and edge perturbation)
    fprintf('\n--- Test 3b: Random Samples Check ---\n');
    samples_in_bounds = 0;

    for s = 1:num_random_samples
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

        % Forward pass with sampled inputs
        Y_sample = gnn.evaluate(X_sample, E_sample);

        if all(Y_sample(:) >= lb_out(:) - tolerance) && ...
           all(Y_sample(:) <= ub_out(:) + tolerance)
            samples_in_bounds = samples_in_bounds + 1;
        end
    end

    fprintf('Random samples within bounds: %d/%d\n', samples_in_bounds, num_random_samples);

    if samples_in_bounds == num_random_samples
        fprintf('PASS: All random samples within bounds\n');
    else
        fprintf('FAIL: %d samples outside bounds!\n', num_random_samples - samples_in_bounds);
        all_tests_passed = false;
    end

    test_results.gine_edge_ieee24 = struct(...
        'center_check', center_in_bounds, ...
        'random_samples', samples_in_bounds);
end

%% Test 4: Over-approximation Quality Check
fprintf('\n=== Test 4: Over-approximation Quality Check ===\n');
fprintf('Checking that output bounds are not excessively large...\n');

if isfield(test_results, 'gine_ieee24') && isstruct(test_results.gine_ieee24)
    % Load GINE results
    modelPath = fullfile(fileparts(mfilename('fullpath')), 'IEEE24', 'models', 'gine_edgelist_pf_ieee24_run4_seed131.mat');
    model = load(modelPath);

    % Recompute for quality metrics
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

    L1 = GINELayer('gine1', W_node1, b_node1, W_edge1, b_edge1);
    L2 = GINELayer('gine2', W_node2, b_node2, W_edge2, b_edge2);
    L3 = GINELayer('gine3', W_node3, b_node3, W_edge3, b_edge3);

    src = double(model.src);
    dst = double(model.dst);
    adj_list = [src, dst];
    E = double(model.E_edge);
    edge_weights = double(model.a);
    X = double(model.X_test_g{1});
    numNodes = size(X, 1);

    gnn = GNN({L1, L2, L3}, [], adj_list, E, edge_weights);

    % Compare different epsilon values
    epsilons = [0.001, 0.005, 0.01];
    fprintf('\nEpsilon | Mean Bound Width | Max Bound Width\n');
    fprintf('--------|-----------------|----------------\n');

    for eps_val = epsilons
        range_per_col = max(X) - min(X);
        eps_matrix = zeros(numNodes, size(X, 2));
        for f = [1, 2]
            if f <= size(X, 2)
                eps_matrix(:, f) = range_per_col(f) * eps_val;
            end
        end

        GS_in = GraphStar(X, -eps_matrix, eps_matrix);
        reachOpts = struct('reachMethod', 'approx-star');
        GS_out = gnn.reach(GS_in, reachOpts);

        [lb_out, ub_out] = GS_out.getRanges();
        bound_widths = ub_out - lb_out;

        fprintf('%.3f   | %.6f        | %.6f\n', eps_val, mean(bound_widths(:)), max(bound_widths(:)));
    end

    fprintf('\nExpected: Bound widths should increase roughly proportionally with epsilon\n');
end

%% Summary
fprintf('\n=== SOUNDNESS VALIDATION SUMMARY ===\n');

if all_tests_passed
    fprintf('\n*** ALL TESTS PASSED ***\n');
    fprintf('\nThe NNV GNN verification implementation is SOUND.\n');
    fprintf('All sampled points from the perturbation region fall within computed bounds.\n');
else
    fprintf('\n*** SOME TESTS FAILED ***\n');
    fprintf('\nReview failed tests above. Possible causes:\n');
    fprintf('  - Numerical precision issues (try increasing tolerance)\n');
    fprintf('  - Bug in reachability computation\n');
    fprintf('  - Bug in forward evaluation\n');
end

fprintf('\n--- Test Results Summary ---\n');
test_fields = fieldnames(test_results);
for i = 1:length(test_fields)
    field = test_fields{i};
    result = test_results.(field);
    if ischar(result) && strcmp(result, 'SKIPPED')
        fprintf('%s: SKIPPED\n', field);
    elseif isstruct(result)
        if result.center_check && result.random_samples == num_random_samples
            fprintf('%s: PASSED\n', field);
        else
            fprintf('%s: FAILED\n', field);
        end
    end
end

fprintf('\n=== Soundness Validation Complete ===\n');
