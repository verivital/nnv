function results = run_gnn_experiments(varargin)
% run_gnn_experiments - GNN verification experiments for CAV26 NNV3.0
%
% This script demonstrates GNN verification capabilities on the IEEE 24-bus
% Power Flow prediction task using three model architectures:
%   - GCN: Graph Convolutional Network
%   - GINE: Graph Isomorphism Network with Edge features (node perturbation)
%   - GINE+Edge: GINE with both node and edge perturbations
%
% Experiments sweep epsilon values [0.001, 0.005, 0.01] and generate
% publication-ready figures for the CAV26 repeatability package.
%
% Usage:
%   run_gnn_experiments()              % Run all experiments, generate figures
%   run_gnn_experiments('no_figures')  % Skip figure generation
%   results = run_gnn_experiments()    % Return results struct
%
% Outputs:
%   results - Struct containing all experiment data and verification results
%
% Author: Anne Tumlin
% Date: 01/15/2026

%% Parse arguments
generate_figures = true;
if nargin > 0 && strcmpi(varargin{1}, 'no_figures')
    generate_figures = false;
end

%% Configuration
models = {'GCN', 'GINE', 'GINE+Edge'};
epsilons = [0.001, 0.005, 0.01];
perturb_features = [1, 2];  % Power injections only (matches GNNV)
v_min = 0.95;  % Voltage spec lower bound (p.u.)
v_max = 1.05;  % Voltage spec upper bound (p.u.)

% Paths (self-contained - all files in this folder)
scriptDir = fileparts(mfilename('fullpath'));
modelDir = fullfile(scriptDir, 'models');
resultsDir = fullfile(scriptDir, 'results');
figuresDir = fullfile(scriptDir, 'figures');

%% Initialize results structure
results = struct();
results.config = struct('models', {models}, 'epsilons', epsilons, ...
    'perturb_features', perturb_features, 'v_min', v_min, 'v_max', v_max);
results.data = cell(length(models), length(epsilons));

%% Print header
fprintf('\n');
fprintf('================================================================\n');
fprintf('        CAV26 GNN Verification Experiments (NNV 3.0)\n');
fprintf('================================================================\n');
fprintf('System: IEEE 24-bus Power Flow\n');
fprintf('Models: %s\n', strjoin(models, ', '));
fprintf('Node epsilon: %s\n', mat2str(epsilons));
fprintf('Edge epsilon: 0.001 (fixed, GINE+Edge only)\n');
fprintf('Voltage spec: [%.2f, %.2f] p.u.\n', v_min, v_max);
fprintf('================================================================\n\n');

total_start = tic;

%% Load models once
fprintf('Loading models...\n');

% Load GINE model
gine_path = fullfile(modelDir, 'gine_ieee24.mat');
gine_model = load(gine_path);

% Load GCN model
gcn_path = fullfile(modelDir, 'gcn_ieee24.mat');
gcn_model = load(gcn_path);

fprintf('Models loaded successfully.\n\n');

%% Run experiments
for m = 1:length(models)
    model_name = models{m};

    for e = 1:length(epsilons)
        epsilon = epsilons(e);

        % Print experiment header (different format for GINE+Edge)
        if strcmp(model_name, 'GINE+Edge')
            fprintf('--- %s, node_eps=%.3f, edge_eps=0.001 ---\n', model_name, epsilon);
        else
            fprintf('--- %s, epsilon=%.3f ---\n', model_name, epsilon);
        end
        exp_start = tic;

        % Run experiment based on model type
        switch model_name
            case 'GCN'
                exp_result = run_gcn_experiment(gcn_model, epsilon, perturb_features, v_min, v_max);
            case 'GINE'
                exp_result = run_gine_experiment(gine_model, epsilon, perturb_features, v_min, v_max);
            case 'GINE+Edge'
                exp_result = run_gine_edge_experiment(gine_model, epsilon, perturb_features, v_min, v_max);
        end

        exp_result.time = toc(exp_start);
        results.data{m, e} = exp_result;

        % Print summary with unknown breakdown
        if exp_result.unknown_boundary > 0 || exp_result.unknown_timeout > 0
            fprintf('  Verified: %d, Violated: %d, Unknown: %d (boundary: %d, timeout: %d) (%.2fs)\n', ...
                exp_result.verified, exp_result.violated, exp_result.unknown, ...
                exp_result.unknown_boundary, exp_result.unknown_timeout, exp_result.time);
        else
            fprintf('  Verified: %d, Violated: %d, Unknown: %d (%.2fs)\n', ...
                exp_result.verified, exp_result.violated, exp_result.unknown, exp_result.time);
        end
    end
    fprintf('\n');
end

total_time = toc(total_start);
fprintf('================================================================\n');
fprintf('Total time: %.2f seconds\n', total_time);
fprintf('================================================================\n\n');

%% Save results
results.total_time = total_time;
results_file = fullfile(resultsDir, 'gnn_results.mat');
save(results_file, 'results');
fprintf('Results saved to: %s\n', results_file);

%% Generate figures
if generate_figures
    fprintf('\nGenerating figures...\n');

    % Original figures (bar chart and line plot)
    generate_model_comparison_figure(results, figuresDir);
    generate_epsilon_sensitivity_figure(results, figuresDir);

    % New domain-specific figures (topology, bounds, dashboard)
    fprintf('Generating domain-specific figures...\n');
    generate_cav26_figures(results, gine_model, figuresDir, 'layout', 'force');

    fprintf('Figures saved to: %s\n', figuresDir);
end

%% Print summary table
fprintf('\n');
fprintf('================================================================\n');
fprintf('                    RESULTS SUMMARY\n');
fprintf('================================================================\n');
fprintf('%-12s | ', 'Model');
for e = 1:length(epsilons)
    fprintf('eps=%.3f | ', epsilons(e));
end
fprintf('\n');
fprintf('%-12s-+-', repmat('-', 1, 12));
for e = 1:length(epsilons)
    fprintf('---------+-');
end
fprintf('\n');

for m = 1:length(models)
    fprintf('%-12s | ', models{m});
    for e = 1:length(epsilons)
        fprintf('  %2d/%2d  | ', results.data{m,e}.verified, ...
            results.data{m,e}.verified + results.data{m,e}.unknown);
    end
    fprintf('\n');
end
fprintf('================================================================\n');

end


%% =========================================================================
%  EXPERIMENT FUNCTIONS
%  =========================================================================

function result = run_gcn_experiment(model, epsilon, perturb_features, v_min, v_max)
% Run GCN verification experiment

    % Extract weights
    params = model.best_params;
    W1 = double(gather(params.mult1.Weights));
    W2 = double(gather(params.mult2.Weights));
    W3 = double(gather(params.mult3.Weights));

    b1 = zeros(size(W1, 2), 1);
    b2 = zeros(size(W2, 2), 1);
    b3 = zeros(size(W3, 2), 1);

    % Create layers with ReLU activations
    L1 = GCNLayer('gcn1', W1, b1);
    R1 = ReluLayer();
    L2 = GCNLayer('gcn2', W2, b2);
    R2 = ReluLayer();
    L3 = GCNLayer('gcn3', W3, b3);
    R3 = ReluLayer();

    % Graph structure
    A_norm = double(model.ANorm_g);
    X = double(model.X_test_g{1});
    numNodes = size(X, 1);

    % Create GNN
    gnn = GNN({L1, R1, L2, R2, L3, R3}, A_norm);

    % Create input perturbation
    GS_in = create_node_perturbation(X, epsilon, perturb_features);

    % Compute reachability (quiet mode)
    reachOpts = struct('reachMethod', 'approx-star');
    t_reach = tic;
    GS_out = gnn.reach(GS_in, reachOpts);
    reach_time = toc(t_reach);

    % Get bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Verify voltage spec
    verif_results = verify_voltage_spec(GS_out, model, v_min, v_max);

    % Extract voltage-specific bounds for figures
    voltage_idx = 3;  % Voltage magnitude index in output
    voltage_lb = lb_out(:, voltage_idx);
    voltage_ub = ub_out(:, voltage_idx);

    % Convert to physical units for visualization
    voltage_lb_phys = voltage_lb * model.global_std_labels(voltage_idx) + model.global_mean_labels(voltage_idx);
    voltage_ub_phys = voltage_ub * model.global_std_labels(voltage_idx) + model.global_mean_labels(voltage_idx);

    % Package results
    result = struct();
    result.reach_time = reach_time;
    result.verified = sum(verif_results == 1);
    result.violated = sum(verif_results == 0);
    result.unknown_boundary = sum(verif_results == 2);
    result.unknown_timeout = sum(verif_results == 3);
    result.unknown = result.unknown_boundary + result.unknown_timeout;
    result.mean_width = mean(ub_out(:) - lb_out(:));
    result.max_width = max(ub_out(:) - lb_out(:));

    % Store per-node data for figures
    result.verif_per_node = verif_results;
    result.voltage_bounds = [voltage_lb_phys, voltage_ub_phys];
    result.bound_widths = ub_out(:, voltage_idx) - lb_out(:, voltage_idx);
end


function result = run_gine_experiment(model, epsilon, perturb_features, v_min, v_max)
% Run GINE verification experiment (node perturbation only)

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

    % Graph structure
    src = double(model.src);
    dst = double(model.dst);
    adj_list = [src, dst];
    E = double(model.E_edge);
    edge_weights = double(model.a);
    X = double(model.X_test_g{1});

    % Create GNN
    gnn = GNN({L1, L2, L3}, [], adj_list, E, edge_weights);

    % Create input perturbation
    GS_in = create_node_perturbation(X, epsilon, perturb_features);

    % Compute reachability (quiet mode)
    reachOpts = struct('reachMethod', 'approx-star');
    t_reach = tic;
    GS_out = gnn.reach(GS_in, reachOpts);
    reach_time = toc(t_reach);

    % Get bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Verify voltage spec
    verif_results = verify_voltage_spec(GS_out, model, v_min, v_max);

    % Extract voltage-specific bounds for figures
    voltage_idx = 3;  % Voltage magnitude index in output
    voltage_lb = lb_out(:, voltage_idx);
    voltage_ub = ub_out(:, voltage_idx);

    % Convert to physical units for visualization
    voltage_lb_phys = voltage_lb * model.global_std_labels(voltage_idx) + model.global_mean_labels(voltage_idx);
    voltage_ub_phys = voltage_ub * model.global_std_labels(voltage_idx) + model.global_mean_labels(voltage_idx);

    % Package results
    result = struct();
    result.reach_time = reach_time;
    result.verified = sum(verif_results == 1);
    result.violated = sum(verif_results == 0);
    result.unknown_boundary = sum(verif_results == 2);
    result.unknown_timeout = sum(verif_results == 3);
    result.unknown = result.unknown_boundary + result.unknown_timeout;
    result.mean_width = mean(ub_out(:) - lb_out(:));
    result.max_width = max(ub_out(:) - lb_out(:));

    % Store per-node data for figures
    result.verif_per_node = verif_results;
    result.voltage_bounds = [voltage_lb_phys, voltage_ub_phys];
    result.bound_widths = ub_out(:, voltage_idx) - lb_out(:, voltage_idx);
end


function result = run_gine_edge_experiment(model, epsilon, perturb_features, v_min, v_max)
% Run GINE verification experiment (node + edge perturbation)

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

    % Graph structure
    src = double(model.src);
    dst = double(model.dst);
    adj_list = [src, dst];
    E = double(model.E_edge);
    edge_weights = double(model.a);
    X = double(model.X_test_g{1});
    numEdges = size(adj_list, 1);

    % Create node perturbation
    GS_in = create_node_perturbation(X, epsilon, perturb_features);

    % Create edge perturbation (fixed at 0.001, independent of node epsilon)
    epsilon_edge = 0.001;  % Fixed edge perturbation
    perturb_edge_features = [1];  % First edge feature (impedance)
    range_per_edge_col = max(E) - min(E);
    eps_matrix_edge = zeros(numEdges, size(E, 2));
    for f = perturb_edge_features
        if f <= size(E, 2)
            eps_matrix_edge(:, f) = range_per_edge_col(f) * epsilon_edge;
        end
    end
    E_star = GraphStar(E, -eps_matrix_edge, eps_matrix_edge);

    % Create GNN with edge perturbation
    gnn = GNN({L1, L2, L3}, [], adj_list, E_star, edge_weights);

    % Compute reachability (quiet mode)
    reachOpts = struct('reachMethod', 'approx-star');
    t_reach = tic;
    GS_out = gnn.reach(GS_in, reachOpts);
    reach_time = toc(t_reach);

    % Get bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Verify voltage spec
    verif_results = verify_voltage_spec(GS_out, model, v_min, v_max);

    % Extract voltage-specific bounds for figures
    voltage_idx = 3;  % Voltage magnitude index in output
    voltage_lb = lb_out(:, voltage_idx);
    voltage_ub = ub_out(:, voltage_idx);

    % Convert to physical units for visualization
    voltage_lb_phys = voltage_lb * model.global_std_labels(voltage_idx) + model.global_mean_labels(voltage_idx);
    voltage_ub_phys = voltage_ub * model.global_std_labels(voltage_idx) + model.global_mean_labels(voltage_idx);

    % Package results
    result = struct();
    result.reach_time = reach_time;
    result.verified = sum(verif_results == 1);
    result.violated = sum(verif_results == 0);
    result.unknown_boundary = sum(verif_results == 2);
    result.unknown_timeout = sum(verif_results == 3);
    result.unknown = result.unknown_boundary + result.unknown_timeout;
    result.mean_width = mean(ub_out(:) - lb_out(:));
    result.max_width = max(ub_out(:) - lb_out(:));

    % Store per-node data for figures
    result.verif_per_node = verif_results;
    result.voltage_bounds = [voltage_lb_phys, voltage_ub_phys];
    result.bound_widths = ub_out(:, voltage_idx) - lb_out(:, voltage_idx);
end


%% =========================================================================
%  HELPER FUNCTIONS
%  =========================================================================

function GS = create_node_perturbation(X, epsilon, perturb_features)
% Create GraphStar with selective feature perturbation

    numNodes = size(X, 1);
    range_per_col = max(X) - min(X);
    eps_matrix = zeros(numNodes, size(X, 2));

    for f = perturb_features
        if f <= size(X, 2)
            eps_matrix(:, f) = range_per_col(f) * epsilon;
        end
    end

    GS = GraphStar(X, -eps_matrix, eps_matrix);
end


%% =========================================================================
%  FIGURE SAVE HELPER
%  =========================================================================

function save_figure(fig, filepath_base)
% Save figure to PDF and PNG without page size warnings

    % Save PNG
    saveas(fig, [filepath_base, '.png']);

    % For PDF, set paper size to match figure
    fig.Units = 'inches';
    fig_pos = fig.Position;
    fig.PaperUnits = 'inches';
    fig.PaperSize = [fig_pos(3), fig_pos(4)];
    fig.PaperPosition = [0, 0, fig_pos(3), fig_pos(4)];

    % Save PDF
    print(fig, [filepath_base, '.pdf'], '-dpdf', '-vector');
end


%% =========================================================================
%  FIGURE GENERATION
%  =========================================================================

function generate_model_comparison_figure(results, figuresDir)
% Generate grouped bar chart comparing models at epsilon=0.01

    models = results.config.models;
    epsilons = results.config.epsilons;

    % Use epsilon=0.01 (index 3)
    eps_idx = find(epsilons == 0.01, 1);
    if isempty(eps_idx)
        eps_idx = length(epsilons);
    end

    % Extract data
    verified = zeros(1, length(models));
    violated = zeros(1, length(models));
    unknown = zeros(1, length(models));

    for m = 1:length(models)
        verified(m) = results.data{m, eps_idx}.verified;
        violated(m) = results.data{m, eps_idx}.violated;
        unknown(m) = results.data{m, eps_idx}.unknown;
    end

    % Create figure
    fig = figure('Position', [100, 100, 600, 400], 'Visible', 'off');

    % Stacked bar chart (same style as dashboard)
    bar_data = [verified; unknown; violated]';
    b = bar(bar_data, 'stacked');

    % Colors (matching dashboard)
    b(1).FaceColor = [0.2, 0.65, 0.3];   % Green for verified
    b(2).FaceColor = [0.95, 0.6, 0.1];   % Orange for unknown
    b(3).FaceColor = [0.8, 0.15, 0.15];  % Red for violated

    % Labels
    set(gca, 'XTickLabel', models);
    xlabel('Model', 'FontSize', 12);
    ylabel('Number of Nodes', 'FontSize', 12);
    title(sprintf('Voltage Specification Verification (\\epsilon=%.3f)', epsilons(eps_idx)), 'FontSize', 14);
    legend({'Verified', 'Unknown', 'Violated'}, 'Location', 'northeast');

    % Grid
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);

    % Save
    save_figure(fig, fullfile(figuresDir, 'gnn_model_comparison'));
    close(fig);
end


function generate_epsilon_sensitivity_figure(results, figuresDir)
% Generate line plot showing verification vs epsilon

    models = results.config.models;
    epsilons = results.config.epsilons;

    % Extract verified counts
    verified_data = zeros(length(models), length(epsilons));
    for m = 1:length(models)
        for e = 1:length(epsilons)
            verified_data(m, e) = results.data{m, e}.verified;
        end
    end

    % Create figure
    fig = figure('Position', [100, 100, 600, 400], 'Visible', 'off');

    colors = lines(length(models));
    markers = {'o-', 's-', 'd-'};

    hold on;
    for m = 1:length(models)
        plot(epsilons, verified_data(m, :), markers{m}, ...
            'Color', colors(m, :), 'LineWidth', 2, 'MarkerSize', 8, ...
            'MarkerFaceColor', colors(m, :));
    end
    hold off;

    % Labels
    xlabel('Perturbation Magnitude (\epsilon)', 'FontSize', 12);
    ylabel('Verified Safe Nodes', 'FontSize', 12);
    title('Verification Sensitivity to Input Perturbation', 'FontSize', 14);
    legend(models, 'Location', 'southwest');

    % Axis settings
    xlim([0, max(epsilons) * 1.1]);
    ylim([0, 15]);

    % Grid
    grid on;
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);

    % Save
    save_figure(fig, fullfile(figuresDir, 'gnn_epsilon_sensitivity'));
    close(fig);
end


%% =========================================================================
%  VERIFICATION FUNCTION (embedded for self-containment)
%  =========================================================================

function results = verify_voltage_spec(GS_out, model_data, v_min, v_max)
% verify_voltage_spec - Verify voltage magnitude bounds on GNN output
%
% Uses LP-based verification (verify_specification) for precise results,
% with fallback to interval bounds when LP returns unknown.

    voltage_idx = 3;   % Index of voltage magnitude in output features
    bus_type_idx = 4;  % Index of bus_type in input features

    % Normalize physical voltage bounds to model space
    v_min_norm = (v_min - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);
    v_max_norm = (v_max - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);

    % Identify voltage-output nodes (bus_type == 1)
    X_physical = model_data.X_test_g{1} .* model_data.global_std + model_data.global_mean;
    voltage_mask = (X_physical(:, bus_type_idx) == 1);

    numNodes = size(GS_out.V, 1);
    numFeatures = size(GS_out.V, 2);
    results = zeros(numNodes, 1);

    % Convert GraphStar to Star for verification
    Y_star = GS_out.toStar();

    for i = 1:numNodes
        if ~voltage_mask(i)
            results(i) = -1;  % Not applicable (not a voltage-output node)
            continue;
        end

        % Extract the i-th node's features as a Star set
        matIdx = zeros(1, numNodes * numFeatures);
        flat_idx = (voltage_idx - 1) * numNodes + i;
        matIdx(flat_idx) = 1;

        % Extract 1D Star for this node's voltage feature
        Y_node = Y_star.affineMap(matIdx, []);

        % Create halfspace constraints for voltage bounds
        G = [1; -1];
        g = [v_max_norm; -v_min_norm];
        Hs = [HalfSpace(G(1,:), g(1)); HalfSpace(G(2,:), g(2))];

        % LP-based verification
        res = verify_specification(Y_node, Hs);

        % Fallback to interval bounds if LP returns unknown
        if res == 2
            [lb, ub] = Y_node.getRanges;
            if lb(1) >= v_min_norm && ub(1) <= v_max_norm
                res = 1;  % Verified safe - bounds fully within spec
            elseif ub(1) < v_min_norm || lb(1) > v_max_norm
                res = 0;  % Violated - bounds fully outside spec
            else
                res = 2;  % Unknown - bounds cross spec boundary
            end
        end

        results(i) = res;
    end
end
