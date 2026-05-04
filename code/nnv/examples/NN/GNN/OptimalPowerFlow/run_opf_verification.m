function results = run_opf_verification(bus_system, layer_type, epsilon, options)
% run_opf_verification - Shared helper for Optimal Power Flow GNN verification
%
% Syntax:
%   results = run_opf_verification(bus_system, layer_type, epsilon)
%   results = run_opf_verification(bus_system, layer_type, epsilon, options)
%
% Inputs:
%   @bus_system: 'IEEE24', 'IEEE39', or 'IEEE118'
%   @layer_type: 'GCN', 'GINE', or 'GINE_edge' (edge perturbation)
%   @epsilon: Perturbation magnitude (e.g., 0.01 for 1%)
%   @options: (optional) struct with fields:
%       perturb_features - Features to perturb [1,2] default (power injections)
%       epsilon_edge - Edge perturbation epsilon (for GINE_edge mode)
%       perturb_edge_features - Edge features to perturb [1] default
%       verbose - Display output (default: true)
%
% Outputs:
%   @results: struct with fields:
%       time - Computation time in seconds
%       GS_out - Output GraphStar
%       lb, ub - Output bounds
%       samples_in_bounds - Number of samples within bounds (out of 10)
%       verified_safe - Number of voltage nodes verified safe
%       unknown - Number of unknown nodes (total)
%       unknown_boundary - Unknown due to bounds crossing spec boundary
%       unknown_timeout - Unknown due to timeout (always 0 in this function)
%       violated - Number of violated nodes
%
% Author: Anne Tumlin
% Date: 01/16/2026

    %% Parse options
    if nargin < 4
        options = struct();
    end

    % Default options
    if ~isfield(options, 'perturb_features')
        options.perturb_features = [1, 2];  % Power injections
    end
    if ~isfield(options, 'epsilon_edge')
        options.epsilon_edge = epsilon;  % Same as node epsilon
    end
    if ~isfield(options, 'perturb_edge_features')
        options.perturb_edge_features = 1;  % Impedance
    end
    if ~isfield(options, 'verbose')
        options.verbose = true;
    end

    verbose = options.verbose;

    %% Validate inputs
    valid_systems = {'IEEE24', 'IEEE39', 'IEEE118'};
    if ~ismember(bus_system, valid_systems)
        error('Invalid bus_system. Must be one of: %s', strjoin(valid_systems, ', '));
    end

    valid_types = {'GCN', 'GINE', 'GINE_edge'};
    if ~ismember(layer_type, valid_types)
        error('Invalid layer_type. Must be one of: %s', strjoin(valid_types, ', '));
    end

    %% Determine model path
    script_dir = fileparts(mfilename('fullpath'));

    % OPF model naming convention (simplified)
    switch layer_type
        case 'GCN'
            model_file = sprintf('gcn_opf_%s.mat', lower(bus_system));
        case {'GINE', 'GINE_edge'}
            model_file = sprintf('gine_opf_%s.mat', lower(bus_system));
    end

    model_path = fullfile(script_dir, bus_system, 'models', model_file);

    if ~exist(model_path, 'file')
        error('Model file not found: %s', model_path);
    end

    if verbose
        fprintf('=== %s Optimal Power Flow Verification (%s) ===\n', layer_type, bus_system);
        fprintf('Model: %s\n', model_path);
        fprintf('Epsilon: %.4f (%.2f%%)\n', epsilon, epsilon*100);
    end

    %% Load model
    model = load(model_path);

    %% Extract graph structure and data
    X = double(model.X_test_g{1});
    numNodes = size(X, 1);

    if strcmp(layer_type, 'GCN')
        % GCN uses normalized adjacency matrix (field name: ANorm_g)
        A_norm = double(model.ANorm_g);
        adj_list = [];
        E = [];
        edge_weights = [];
    else
        % GINE uses edge list
        src = double(model.src);
        dst = double(model.dst);
        adj_list = [src, dst];
        E = double(model.E_edge);
        edge_weights = double(model.a);
        A_norm = [];
    end

    numEdges = size(adj_list, 1);

    if verbose
        fprintf('Graph: %d nodes', numNodes);
        if numEdges > 0
            fprintf(', %d edges', numEdges);
        end
        fprintf('\n');
    end

    %% Create layers
    params = model.best_params;

    if strcmp(layer_type, 'GCN')
        % GCN layers
        W1 = double(gather(params.mult1.Weights));
        W2 = double(gather(params.mult2.Weights));
        W3 = double(gather(params.mult3.Weights));
        b1 = zeros(size(W1, 2), 1);
        b2 = zeros(size(W2, 2), 1);
        b3 = zeros(size(W3, 2), 1);

        L1 = GCNLayer('gcn1', W1, b1);
        L2 = GCNLayer('gcn2', W2, b2);
        L3 = GCNLayer('gcn3', W3, b3);

        gnn = GNN({L1, L2, L3}, A_norm);
    else
        % GINE layers
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

        % Handle edge perturbation mode
        if strcmp(layer_type, 'GINE_edge')
            % Create edge perturbation GraphStar
            range_edge = max(E) - min(E);
            eps_edge = zeros(size(E));
            for f = options.perturb_edge_features
                if f <= size(E, 2)
                    eps_edge(:, f) = range_edge(f) * options.epsilon_edge;
                end
            end
            E_star = GraphStar(E, -eps_edge, eps_edge);
            gnn = GNN({L1, L2, L3}, [], adj_list, E_star, edge_weights);
        else
            gnn = GNN({L1, L2, L3}, [], adj_list, E, edge_weights);
        end
    end

    %% Create node perturbation
    range_per_col = max(X) - min(X);
    eps_matrix = zeros(numNodes, size(X, 2));

    for f = options.perturb_features
        if f <= size(X, 2)
            eps_matrix(:, f) = range_per_col(f) * epsilon;
        end
    end

    GS_in = GraphStar(X, -eps_matrix, eps_matrix);

    %% Compute reachability
    if verbose
        fprintf('\nComputing reachability...\n');
    end

    reachOpts = struct('reachMethod', 'approx-star');
    if verbose
        reachOpts.dis_opt = 'display';
    end

    t_start = tic;
    GS_out = gnn.reach(GS_in, reachOpts);
    comp_time = toc(t_start);

    if verbose
        fprintf('Completed in %.2f seconds\n', comp_time);
    end

    %% Get output bounds
    [lb, ub] = GS_out.getRanges();

    %% Sample validation
    num_samples = 10;
    samples_in_bounds = 0;

    for s = 1:num_samples
        alpha = rand(GS_in.numPred, 1) .* (GS_in.pred_ub - GS_in.pred_lb) + GS_in.pred_lb;
        X_sample = GS_in.V(:, :, 1);
        for k = 1:GS_in.numPred
            X_sample = X_sample + alpha(k) * GS_in.V(:, :, k+1);
        end

        if strcmp(layer_type, 'GINE_edge')
            % Sample edge features too
            alpha_edge = rand(E_star.numPred, 1) .* (E_star.pred_ub - E_star.pred_lb) + E_star.pred_lb;
            E_sample = E_star.V(:, :, 1);
            for k = 1:E_star.numPred
                E_sample = E_sample + alpha_edge(k) * E_star.V(:, :, k+1);
            end
            Y_sample = gnn.evaluate(X_sample, E_sample);
        else
            Y_sample = gnn.evaluate(X_sample);
        end

        tol = 1e-6;
        if all(Y_sample(:) >= lb(:) - tol) && all(Y_sample(:) <= ub(:) + tol)
            samples_in_bounds = samples_in_bounds + 1;
        end
    end

    if verbose
        fprintf('Samples within bounds: %d/%d\n', samples_in_bounds, num_samples);
    end

    %% Voltage verification
    [verified_safe, unknown_boundary, violated] = verify_voltage_nodes(GS_out, model, verbose);

    %% Package results
    results = struct();
    results.time = comp_time;
    results.GS_out = GS_out;
    results.lb = lb;
    results.ub = ub;
    results.samples_in_bounds = samples_in_bounds;
    results.verified_safe = verified_safe;
    results.unknown = unknown_boundary;  % Total unknown (for backwards compatibility)
    results.unknown_boundary = unknown_boundary;  % Bounds cross spec boundary
    results.unknown_timeout = 0;  % Timeout not applicable in this function
    results.violated = violated;
    results.bus_system = bus_system;
    results.layer_type = layer_type;
    results.epsilon = epsilon;
end


function [verified_safe, unknown_boundary, violated] = verify_voltage_nodes(GS_out, model, verbose)
% Verify voltage specification for all voltage nodes
% Returns counts of: verified_safe, unknown_boundary (bounds cross spec), violated

    % Voltage specification (per-unit)
    v_min = 0.95;
    v_max = 1.05;

    % Get normalization parameters
    global_mean = model.global_mean_labels;
    global_std = model.global_std_labels;

    % Voltage is feature index 3
    voltage_idx = 3;
    v_mean = global_mean(voltage_idx);
    v_std = global_std(voltage_idx);

    % Normalized spec bounds
    spec_lb = (v_min - v_mean) / v_std;
    spec_ub = (v_max - v_mean) / v_std;

    % Get output bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Get input data for bus type
    X = double(model.X_test_g{1});
    X_physical = X .* model.global_std + model.global_mean;
    bus_type = round(X_physical(:, 4));

    % Count by status
    verified_safe = 0;
    unknown_boundary = 0;  % Bounds cross spec boundary
    violated = 0;

    numNodes = size(lb_out, 1);

    for n = 1:numNodes
        if bus_type(n) ~= 1
            continue;  % Skip non-voltage buses
        end

        node_lb = lb_out(n, voltage_idx);
        node_ub = ub_out(n, voltage_idx);

        if node_lb >= spec_lb && node_ub <= spec_ub
            verified_safe = verified_safe + 1;
        elseif node_ub < spec_lb || node_lb > spec_ub
            violated = violated + 1;
        else
            unknown_boundary = unknown_boundary + 1;  % Bounds cross spec boundary
        end
    end

    if verbose
        fprintf('\n=== Voltage Verification ===\n');
        fprintf('Spec: %.2f <= V <= %.2f p.u.\n', v_min, v_max);
        fprintf('  Verified safe: %d nodes\n', verified_safe);
        fprintf('  Violated: %d nodes\n', violated);
        fprintf('  Unknown: %d nodes\n', unknown_boundary);
        if unknown_boundary > 0
            fprintf('    - Bounds cross spec boundary: %d\n', unknown_boundary);
        end
    end
end
