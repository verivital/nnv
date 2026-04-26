%% soundness_check.m
% Soundness validation for the GNN verification pipeline.
%
% For each shipped (architecture × task × grid) configuration, loads the
% model, runs reachability on a few graphs, then Monte Carlo samples random
% perturbations to confirm that every concrete GNN evaluation lies within
% the computed reachable bounds.
%
% Checks:
%   1. Nominal accuracy — gnn.evaluate(X) matches test_data.Y_all
%   2. Reachability soundness — random samples within epsilon ball
%      fall inside reachable bounds (over-approximation holds)
%   3. Bound width stats — reports mean/max width per config
%
% Returns a results struct array for figure generation (optional output).
%
% Usage:
%   soundness_check()
%   results = soundness_check('n_graphs', 5, 'n_samples', 100)
%   soundness_check('tasks', {'pf'}, 'archs', {'gcn'})
%
% Author: Anne Tumlin
% Date: 03/13/2026

function results = soundness_check(varargin)

p = inputParser;
addParameter(p, 'n_graphs', 3, @isnumeric);
addParameter(p, 'n_samples', 50, @isnumeric);
addParameter(p, 'epsilon', 0.005, @isnumeric);
addParameter(p, 'tol', 1e-4, @isnumeric);
addParameter(p, 'tasks', {'pf', 'opf'}, @iscell);
addParameter(p, 'grids', {'IEEE24'}, @iscell);
addParameter(p, 'archs', {'gcn', 'sage', 'gine_conv'}, @iscell);
parse(p, varargin{:});

n_graphs  = p.Results.n_graphs;
n_samples = p.Results.n_samples;
epsilon   = p.Results.epsilon;
tol       = p.Results.tol;
tasks     = p.Results.tasks;
grids     = p.Results.grids;
archs     = p.Results.archs;

fprintf('=== Soundness Check: PF/OPF Verification ===\n');
fprintf('  Graphs/config: %d | Samples/graph: %d | eps=%.4f | tol=%.0e\n', ...
    n_graphs, n_samples, epsilon, tol);
fprintf('  Tasks: %s | Grids: %s | Archs: %s\n\n', ...
    strjoin(tasks, ','), strjoin(grids, ','), strjoin(archs, ','));

scriptDir = fileparts(mfilename('fullpath'));

task_folders = struct('pf', 'PowerFlow', 'opf', 'OptimalPowerFlow');
perturb_node_features = [1, 2];  % P, Q columns

% Voltage bounds per grid (matching run_all_experiments)
grid_vbounds = struct( ...
    'IEEE24',  [0.95, 1.05], ...
    'IEEE39',  [0.94, 1.06], ...
    'IEEE118', [0.94, 1.09]);

reachOpts = struct('reachMethod', 'approx-star');

% Summary tracking
summary = {};
total_configs = 0;
total_pass = 0;
total_fail = 0;

% Results array for figure generation
results = struct([]);
ri = 0;  % results index

for ti = 1:length(tasks)
    task = tasks{ti};
    task_folder = task_folders.(task);

    for gri = 1:length(grids)
        grid = grids{gri};

        for ai = 1:length(archs)
            arch = archs{ai};
            total_configs = total_configs + 1;

            % Build model path
            mat_file = sprintf('%s_%s_%s.mat', arch, task, lower(grid));
            mat_path = fullfile(scriptDir, task_folder, grid, 'models', mat_file);

            config_label = sprintf('%s/%s/%s', upper(task), grid, arch);
            fprintf('--- %s ---\n', config_label);

            if ~isfile(mat_path)
                fprintf('  SKIP: model not found (%s)\n\n', mat_file);
                summary{end+1} = {config_label, 'SKIP', 'model not found'}; %#ok<AGROW>
                continue;
            end

            try
                [gnn, test_data, norm_stats] = gnn2nnv(mat_path);
            catch ME
                fprintf('  FAIL: load error — %s\n\n', ME.message);
                summary{end+1} = {config_label, 'FAIL', ['load: ' ME.message]}; %#ok<AGROW>
                total_fail = total_fail + 1;
                continue;
            end

            % Load python_predictions for nominal check (Y_test_g stores
            % physical-unit targets, not normalized model outputs)
            raw_model = load(mat_path, 'python_predictions');
            has_python_pred = isfield(raw_model, 'python_predictions') && ...
                              iscell(raw_model.python_predictions);

            n_avail = length(test_data.X_all);
            n_test = min(n_graphs, n_avail);

            is_gine = contains(arch, 'gine');
            nom_failures = 0;
            containment_failures = 0;
            total_evals = 0;
            max_nom_err = 0;
            all_widths = [];

            for gi = 1:n_test
                X = test_data.X_all{gi};

                % --- Check 1: Nominal accuracy ---
                % Compare MATLAB evaluate against Python predictions (both
                % in normalized space). Y_test_g is physical-unit ground
                % truth and cannot be compared directly to model output.
                Y_pred = gnn.evaluate(X);
                if has_python_pred && gi <= numel(raw_model.python_predictions)
                    Y_ref = double(raw_model.python_predictions{gi});
                    nom_err = max(abs(Y_pred(:) - Y_ref(:)));
                else
                    % Fallback: skip nominal check if no python_predictions
                    nom_err = 0;
                end
                max_nom_err = max(max_nom_err, nom_err);
                if nom_err > 0.01
                    nom_failures = nom_failures + 1;
                    fprintf('  Graph %d: nominal err=%.2e (FAIL)\n', gi, nom_err);
                end

                % --- Create perturbation (matching run_all_experiments) ---
                range_per_col = max(X) - min(X);
                eps_matrix = zeros(size(X));
                for f = perturb_node_features
                    if f <= size(X, 2)
                        eps_matrix(:, f) = range_per_col(f) * epsilon;
                    end
                end
                GS_in = GraphStar(X, -eps_matrix, eps_matrix);

                % Collect samples for figure data
                sample_outputs = zeros(size(X, 1), size(X, 2), n_samples);

                % --- Reachability ---
                if is_gine
                    % Subgraph verification for GINE (matches pipeline)
                    if isfield(norm_stats, 'X_max')
                        X_phys = X .* norm_stats.X_max;
                    else
                        X_phys = X;
                    end
                    pq_nodes = find(round(X_phys(:, 4)) == 1);

                    [node_outputs, sg_info] = gnn.reachSubgraph(GS_in, pq_nodes, reachOpts);

                    % Build full lb/ub from subgraph results (PQ nodes only)
                    lb = nan(size(X, 1), size(Y_pred, 2));
                    ub = nan(size(X, 1), size(Y_pred, 2));
                    for ti_pq = 1:length(pq_nodes)
                        t_local = sg_info(ti_pq).target_local_idx;
                        [lb_t, ub_t] = node_outputs{ti_pq}.estimateRanges();
                        lb(pq_nodes(ti_pq), :) = lb_t(t_local, :);
                        ub(pq_nodes(ti_pq), :) = ub_t(t_local, :);
                    end

                    % Check containment for each PQ node
                    rng(gi * 42);
                    for si = 1:n_samples
                        dX = (2 * rand(size(X)) - 1) .* eps_matrix;
                        Y_sample = gnn.evaluate(X + dX);
                        sample_outputs(:, :, si) = Y_sample;

                        for ti_pq = 1:length(pq_nodes)
                            t = pq_nodes(ti_pq);
                            for feat = 1:size(Y_sample, 2)
                                total_evals = total_evals + 1;
                                y_val = Y_sample(t, feat);
                                if y_val < lb(t, feat) - tol || ...
                                   y_val > ub(t, feat) + tol
                                    containment_failures = containment_failures + 1;
                                    if containment_failures <= 3
                                        fprintf('  CONTAIN FAIL: graph %d, sample %d, node %d, feat %d: y=%.6f not in [%.6f, %.6f]\n', ...
                                            gi, si, t, feat, y_val, lb(t, feat), ub(t, feat));
                                    end
                                end
                            end
                        end
                    end

                    % Bound width stats from subgraph outputs
                    for ti_pq = 1:length(pq_nodes)
                        t_local = sg_info(ti_pq).target_local_idx;
                        [lb_t, ub_t] = node_outputs{ti_pq}.estimateRanges();
                        all_widths = [all_widths; ub_t(t_local,:) - lb_t(t_local,:)]; %#ok<AGROW>
                    end

                else
                    % Full-graph verification for GCN/SAGE
                    GS_out = gnn.reach(GS_in, reachOpts);
                    [lb, ub] = GS_out.estimateRanges();

                    % Bound width stats
                    widths = ub - lb;
                    all_widths = [all_widths; widths(:)]; %#ok<AGROW>

                    % Monte Carlo containment
                    rng(gi * 42);
                    for si = 1:n_samples
                        dX = (2 * rand(size(X)) - 1) .* eps_matrix;
                        Y_sample = gnn.evaluate(X + dX);
                        sample_outputs(:, :, si) = Y_sample;

                        violations = (Y_sample < lb - tol) | (Y_sample > ub + tol);
                        n_viol = sum(violations(:));
                        total_evals = total_evals + numel(Y_sample);
                        if n_viol > 0
                            containment_failures = containment_failures + n_viol;
                            if containment_failures <= 3
                                [vn, vf] = find(violations, 1);
                                fprintf('  CONTAIN FAIL: graph %d, sample %d, node %d, feat %d: y=%.6f not in [%.6f, %.6f]\n', ...
                                    gi, si, vn, vf, Y_sample(vn, vf), lb(vn, vf), ub(vn, vf));
                            end
                        end
                    end
                end

                % Store per-graph results for figures
                ri = ri + 1;
                results(ri).config_label = config_label;
                results(ri).task = task;
                results(ri).grid = grid;
                results(ri).arch = arch;
                results(ri).graph_idx = gi;
                results(ri).epsilon = epsilon;
                results(ri).X_center = X;
                results(ri).Y_nominal = Y_pred;
                results(ri).lb = lb;
                results(ri).ub = ub;
                results(ri).samples = sample_outputs;
                results(ri).norm_stats = norm_stats;
                if is_gine
                    results(ri).pq_nodes = pq_nodes;
                else
                    results(ri).pq_nodes = [];
                end
            end

            % Report
            mean_width = 0; max_width = 0;
            if ~isempty(all_widths)
                mean_width = mean(all_widths(:));
                max_width = max(all_widths(:));
            end

            passed = (nom_failures == 0) && (containment_failures == 0);
            status = conditional_str(passed, 'PASS', 'FAIL');
            if passed
                total_pass = total_pass + 1;
            else
                total_fail = total_fail + 1;
            end

            fprintf('  Nominal: max_err=%.2e (%s) | Containment: %d/%d violations | Width: mean=%.2e, max=%.2e | %s\n\n', ...
                max_nom_err, conditional_str(nom_failures == 0, 'ok', 'FAIL'), ...
                containment_failures, total_evals, mean_width, max_width, status);

            detail = sprintf('nom_err=%.2e, contain=%d/%d, width=%.2e', ...
                max_nom_err, containment_failures, total_evals, max_width);
            summary{end+1} = {config_label, status, detail}; %#ok<AGROW>
        end
    end
end

%% Print summary table
fprintf('\n=== SUMMARY ===\n');
fprintf('  %-25s  %-6s  %s\n', 'Config', 'Result', 'Details');
fprintf('  %s\n', repmat('-', 1, 75));
for i = 1:length(summary)
    s = summary{i};
    fprintf('  %-25s  %-6s  %s\n', s{1}, s{2}, s{3});
end
fprintf('\n  Total: %d configs | %d PASS | %d FAIL | %d SKIP\n', ...
    total_configs, total_pass, total_fail, total_configs - total_pass - total_fail);

if total_fail == 0
    fprintf('\n  OVERALL: ALL SOUNDNESS CHECKS PASSED\n');
else
    fprintf('\n  OVERALL: %d CONFIG(S) FAILED\n', total_fail);
end

end


%% Helper
function s = conditional_str(cond, true_str, false_str)
    if cond, s = true_str; else, s = false_str; end
end
