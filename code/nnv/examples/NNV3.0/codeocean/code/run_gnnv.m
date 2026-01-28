function run_gnnv()
%% run_gnnv - GNNV (Graph Neural Network Verification) test for CodeOcean
%
% This runs GNN verification experiments on the IEEE 24-bus Power Flow
% prediction task using GCN and GINE model architectures.
%
% Models: GCN, GINE, GINE+Edge
% Epsilon values: [0.001, 0.005, 0.01]
% Voltage spec: [0.95, 1.05] p.u.

disp('Running GNNV (Graph Neural Network Verification) test...');

%% Configuration
models = {'GCN', 'GINE', 'GINE+Edge'};
epsilons = [0.001, 0.005, 0.01];
perturb_features = [1, 2];  % Power injections only
v_min = 0.95;  % Voltage spec lower bound (p.u.)
v_max = 1.05;  % Voltage spec upper bound (p.u.)

% Use 10 scenarios to match original NNV 3.0 experiments
num_scenarios = 10;
scenario_indices = 1:100:1000;  % Evenly-spaced: [1, 101, 201, ..., 901]

% CodeOcean paths
modelDir = '/data/GNNV/models';
resultsDir = '/results/GNNV';

% Ensure results directory exists
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% Initialize results structure
results = struct();
results.config = struct('models', {models}, 'epsilons', epsilons, ...
    'perturb_features', perturb_features, 'v_min', v_min, 'v_max', v_max, ...
    'num_scenarios', num_scenarios, 'scenario_indices', scenario_indices);
results.data = cell(length(models), length(epsilons), num_scenarios);

%% Print header
fprintf('\n');
fprintf('================================================================\n');
fprintf('        GNNV Verification Experiments (NNV 3.0)\n');
fprintf('================================================================\n');
fprintf('System: IEEE 24-bus Power Flow\n');
fprintf('Models: %s\n', strjoin(models, ', '));
fprintf('Node epsilon: %s\n', mat2str(epsilons));
fprintf('Edge epsilon: 0.001 (fixed, GINE+Edge only)\n');
fprintf('Voltage spec: [%.2f, %.2f] p.u.\n', v_min, v_max);
fprintf('Scenarios: %d (indices: %s)\n', num_scenarios, mat2str(scenario_indices));
fprintf('Total experiments: %d\n', length(models) * length(epsilons) * num_scenarios);
fprintf('================================================================\n\n');

total_start = tic;

%% Load models
fprintf('Loading models...\n');

% Load GINE model
gine_path = fullfile(modelDir, 'gine_ieee24.mat');
if ~exist(gine_path, 'file')
    error('GINE model not found: %s', gine_path);
end
gine_model = load(gine_path);

% Load GCN model
gcn_path = fullfile(modelDir, 'gcn_ieee24.mat');
if ~exist(gcn_path, 'file')
    error('GCN model not found: %s', gcn_path);
end
gcn_model = load(gcn_path);

fprintf('Models loaded successfully.\n\n');

%% Run experiments
exp_count = 0;
total_experiments = length(models) * length(epsilons) * num_scenarios;

for s = 1:num_scenarios
    scenario_idx = scenario_indices(s);
    fprintf('=== Scenario %d/%d (sample index: %d) ===\n', s, num_scenarios, scenario_idx);

    for m = 1:length(models)
        model_name = models{m};

        for e = 1:length(epsilons)
            epsilon = epsilons(e);
            exp_count = exp_count + 1;

            % Print experiment header
            if strcmp(model_name, 'GINE+Edge')
                fprintf('[%d/%d] %s, eps=%.3f, edge_eps=0.001: ', exp_count, total_experiments, model_name, epsilon);
            else
                fprintf('[%d/%d] %s, eps=%.3f: ', exp_count, total_experiments, model_name, epsilon);
            end
            exp_start = tic;

            % Run experiment based on model type
            switch model_name
                case 'GCN'
                    exp_result = run_gcn_experiment(gcn_model, epsilon, perturb_features, v_min, v_max, scenario_idx);
                case 'GINE'
                    exp_result = run_gine_experiment(gine_model, epsilon, perturb_features, v_min, v_max, scenario_idx);
                case 'GINE+Edge'
                    exp_result = run_gine_edge_experiment(gine_model, epsilon, perturb_features, v_min, v_max, scenario_idx);
            end

            exp_result.time = toc(exp_start);
            exp_result.scenario_idx = scenario_idx;
            results.data{m, e, s} = exp_result;

            % Print compact summary
            fprintf('V=%d, X=%d, U=%d (%.2fs)\n', ...
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

%% Print summary table
fprintf('\n');
fprintf('================================================================\n');
fprintf('              RESULTS SUMMARY (mean over %d scenarios)\n', num_scenarios);
fprintf('================================================================\n');
fprintf('%-12s | ', 'Model');
for e = 1:length(epsilons)
    fprintf('  eps=%.3f   | ', epsilons(e));
end
fprintf('\n');
fprintf('%-12s-+-', repmat('-', 1, 12));
for e = 1:length(epsilons)
    fprintf('-------------+-');
end
fprintf('\n');

for m = 1:length(models)
    fprintf('%-12s | ', models{m});
    for e = 1:length(epsilons)
        % Collect verified counts across all scenarios
        verified_counts = zeros(1, num_scenarios);
        for ss = 1:num_scenarios
            verified_counts(ss) = results.data{m, e, ss}.verified;
        end
        mean_v = mean(verified_counts);
        std_v = std(verified_counts);
        fprintf('%5.1f ± %4.1f | ', mean_v, std_v);
    end
    fprintf('\n');
end
fprintf('================================================================\n');

% Print total nodes for reference
fprintf('Note: Total voltage-output nodes per scenario: %d\n', ...
    results.data{1,1,1}.verified + results.data{1,1,1}.unknown + results.data{1,1,1}.violated);

disp('GNNV test completed successfully.');

end


%% =========================================================================
%  EXPERIMENT FUNCTIONS
%  =========================================================================

function result = run_gcn_experiment(model, epsilon, perturb_features, v_min, v_max, scenario_idx)
% Run GCN verification experiment

    % Extract weights
    params = model.best_params;
    W1 = double(gather(params.mult1.Weights));
    W2 = double(gather(params.mult2.Weights));
    W3 = double(gather(params.mult3.Weights));

    % Extract bias from model (if present), otherwise use zeros
    if isfield(params.mult1, 'Bias')
        b1 = double(extractdata(gather(params.mult1.Bias)));
        b2 = double(extractdata(gather(params.mult2.Bias)));
        b3 = double(extractdata(gather(params.mult3.Bias)));
    else
        b1 = zeros(size(W1, 2), 1);
        b2 = zeros(size(W2, 2), 1);
        b3 = zeros(size(W3, 2), 1);
    end

    % Create layers with ReLU activations
    L1 = GCNLayer('gcn1', W1, b1);
    R1 = ReluLayer();
    L2 = GCNLayer('gcn2', W2, b2);
    R2 = ReluLayer();
    L3 = GCNLayer('gcn3', W3, b3);
    R3 = ReluLayer();

    % Graph structure
    A_norm = double(model.ANorm_g);
    X = double(model.X_test_g{scenario_idx});

    % Create GNN
    gnn = GNN({L1, R1, L2, R2, L3, R3}, A_norm);

    % Create input perturbation
    GS_in = create_node_perturbation(X, epsilon, perturb_features);

    % Compute reachability
    reachOpts = struct('reachMethod', 'approx-star');
    t_reach = tic;
    GS_out = gnn.reach(GS_in, reachOpts);
    reach_time = toc(t_reach);

    % Get bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Verify voltage spec
    verif_results = verify_voltage_spec(GS_out, model, v_min, v_max, scenario_idx);

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
    result.verif_per_node = verif_results;
end


function result = run_gine_experiment(model, epsilon, perturb_features, v_min, v_max, scenario_idx)
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
    X = double(model.X_test_g{scenario_idx});

    % Create GNN
    gnn = GNN({L1, L2, L3}, [], adj_list, E, edge_weights);

    % Create input perturbation
    GS_in = create_node_perturbation(X, epsilon, perturb_features);

    % Compute reachability
    reachOpts = struct('reachMethod', 'approx-star');
    t_reach = tic;
    GS_out = gnn.reach(GS_in, reachOpts);
    reach_time = toc(t_reach);

    % Get bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Verify voltage spec
    verif_results = verify_voltage_spec(GS_out, model, v_min, v_max, scenario_idx);

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
    result.verif_per_node = verif_results;
end


function result = run_gine_edge_experiment(model, epsilon, perturb_features, v_min, v_max, scenario_idx)
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
    X = double(model.X_test_g{scenario_idx});
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

    % Compute reachability
    reachOpts = struct('reachMethod', 'approx-star');
    t_reach = tic;
    GS_out = gnn.reach(GS_in, reachOpts);
    reach_time = toc(t_reach);

    % Get bounds
    [lb_out, ub_out] = GS_out.getRanges();

    % Verify voltage spec
    verif_results = verify_voltage_spec(GS_out, model, v_min, v_max, scenario_idx);

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
    result.verif_per_node = verif_results;
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


function results = verify_voltage_spec(GS_out, model_data, v_min, v_max, scenario_idx)
% verify_voltage_spec - Verify voltage magnitude bounds on GNN output

    voltage_idx = 3;   % Index of voltage magnitude in output features
    bus_type_idx = 4;  % Index of bus_type in input features

    % Normalize physical voltage bounds to model space
    v_min_norm = (v_min - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);
    v_max_norm = (v_max - model_data.global_mean_labels(voltage_idx)) / ...
                 model_data.global_std_labels(voltage_idx);

    % Identify voltage-output nodes (bus_type == 1)
    X_physical = model_data.X_test_g{scenario_idx} .* model_data.global_std + model_data.global_mean;
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
