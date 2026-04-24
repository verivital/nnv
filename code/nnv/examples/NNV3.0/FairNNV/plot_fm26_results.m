%% Plot FM26 FairNNV Results
% Generates figures and LaTeX tables from FM26 verification results
% Outputs:
%   (1) LaTeX table for counterfactual fairness
%   (2) Combined individual fairness area plot (smooth lines with filled regions)
%   (3) LaTeX table for timing results (separated by fairness type)
%
% This script can be run standalone or called from run_fm26_fairnnv.m
% When run standalone, default paths are used.
% When called from runner, paths are set via config struct.

%% Setup
% Check if config exists (set by runner script), otherwise use defaults
if ~exist('config', 'var')
    % Default configuration for standalone execution
    config.outputDir = './fm26_fairnnv_results';
    config.savePNG = true;
    config.savePDF = true;
end

resultsDir = config.outputDir;

% Find the most recent CSV files in the results directory
counterfactual_files = dir(fullfile(resultsDir, 'fm26_counterfactual_*.csv'));
individual_files = dir(fullfile(resultsDir, 'fm26_individual_*.csv'));
timing_files = dir(fullfile(resultsDir, 'fm26_timing_*.csv'));

% Check if files exist
if isempty(counterfactual_files) || isempty(individual_files) || isempty(timing_files)
    error('CSV files not found in %s. Please run adult_verify_fm26.m first.', resultsDir);
end

% Get the most recent files (sorted by date)
[~, idx] = max([counterfactual_files.datenum]);
csv_counterfactual = fullfile(resultsDir, counterfactual_files(idx).name);

[~, idx] = max([individual_files.datenum]);
csv_individual = fullfile(resultsDir, individual_files(idx).name);

[~, idx] = max([timing_files.datenum]);
csv_timing = fullfile(resultsDir, timing_files(idx).name);

disp("Loading results from:");
disp("  " + csv_counterfactual);
disp("  " + csv_individual);
disp("  " + csv_timing);

%% Load CSV Data
% Counterfactual fairness data
counterfactual_data = readtable(csv_counterfactual);

% Individual fairness data
individual_data = readtable(csv_individual);

% Timing data
timing_data = readtable(csv_timing);

%% Figure Settings
% Professional color scheme (matching Python matplotlib style)
% Hex colors: #2ecc71 (green), #e74c3c (red)
color_fair = [46/255, 204/255, 113/255];     % Green (#2ecc71)
color_unfair = [231/255, 76/255, 60/255];    % Red (#e74c3c)

% Get unique models
models = unique(individual_data.Model);

% Model display names (fuller titles for figures/tables)
% Based on actual architectures:
% AC-1: 13→16→8→2 (~350 params) - Small
% AC-3: 13→50→2 (~750 params) - Medium
% AC-4: 13→100→100→2 (~11,500 params) - Large
model_display_names = containers.Map();
model_display_names('AC-1') = 'Adult Census - Small Model';
model_display_names('AC-3') = 'Adult Census - Medium Model';
model_display_names('AC-4') = 'Adult Census - Large Model';

%% LaTeX Table 1: Counterfactual Fairness (with timing)
disp(" ");
disp("======= COUNTERFACTUAL FAIRNESS LATEX TABLE ==========");

latex_cf_filename = fullfile(resultsDir, 'fm26_counterfactual_table.tex');
fid = fopen(latex_cf_filename, 'w');

fprintf(fid, '\\begin{table}[ht]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Counterfactual Fairness Verification Results ($\\epsilon = 0$)}\n');
fprintf(fid, '\\label{tab:counterfactual_fairness}\n');
fprintf(fid, '\\begin{tabular}{lccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Model & Fair (\\%%) & Unfair (\\%%) & Avg. Time (s) \\\\\n');
fprintf(fid, '\\midrule\n');

for r = 1:height(counterfactual_data)
    modelName = counterfactual_data.Model{r};
    % Use display name if available
    if isKey(model_display_names, modelName)
        displayName = model_display_names(modelName);
    else
        displayName = modelName;
    end
    % Find corresponding timing (epsilon = 0)
    timing_idx = strcmp(timing_data.Model, modelName) & timing_data.Epsilon == 0;
    if any(timing_idx)
        avg_time = timing_data.AvgTimePerSample(timing_idx);
    else
        avg_time = NaN;
    end
    fprintf(fid, '%s & %.1f & %.1f & %.4f \\\\\n', ...
        displayName, ...
        counterfactual_data.FairPercent(r), ...
        counterfactual_data.UnfairPercent(r), ...
        avg_time);
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);

% Also print to console
type(latex_cf_filename);
disp(" ");
disp("Saved: " + latex_cf_filename);

%% Figure: Combined Individual Fairness (Stacked Area Plot - Professional Style)
n_models = length(models);
fig = figure('Name', 'Individual Fairness - All Models', 'Position', [100, 100, 600*n_models, 500]);

for m = 1:n_models
    ax = subplot(1, n_models, m);

    modelName = models{m};

    % Filter data for this model
    model_data = individual_data(strcmp(individual_data.Model, modelName), :);

    % Sort by epsilon
    model_data = sortrows(model_data, 'Epsilon');

    % Extract values
    epsilons = model_data.Epsilon;
    fair_pct = model_data.FairPercent;
    unfair_pct = model_data.UnfairPercent;
    n_eps = length(epsilons);

    % Use discrete x positions (like Python version)
    x = 0:(n_eps-1);

    % Create stacked area plot (bottom to top: unfair, then fair)
    % This matches the Python stacking order
    y1 = unfair_pct';           % Bottom layer: Unfair
    y2 = y1 + fair_pct';        % Top layer: Fair (should sum to 100)

    hold on;

    % Fill area for Unfair (bottom, from 0 to y1)
    fill([x, fliplr(x)], [y1, zeros(1, n_eps)], ...
         color_unfair, 'FaceAlpha', 0.9, 'EdgeColor', 'none');

    % Fill area for Fair (top, from y1 to y2)
    fill([x, fliplr(x)], [y2, fliplr(y1)], ...
         color_fair, 'FaceAlpha', 0.9, 'EdgeColor', 'none');

    % Add white edge line between areas for clarity
    plot(x, y1, 'w-', 'LineWidth', 1.5);

    hold off;

    % Labels and formatting with bold fonts and larger size
    xlabel('Perturbation Level (\epsilon)', 'FontWeight', 'bold', 'FontSize', 12);
    ylabel('Percentage (%)', 'FontWeight', 'bold', 'FontSize', 12);

    % Use full display name for title
    if isKey(model_display_names, modelName)
        displayTitle = model_display_names(modelName);
    else
        displayTitle = modelName;
    end
    title(displayTitle, 'FontWeight', 'bold', 'FontSize', 14);

    % Set x-axis with epsilon labels
    set(ax, 'XTick', x);
    eps_labels = arrayfun(@(e) sprintf('%.2f', e), epsilons, 'UniformOutput', false);
    set(ax, 'XTickLabel', eps_labels);

    % Axis limits
    ylim([0 100]);
    xlim([-0.3, n_eps - 0.7]);

    % Professional grid styling
    grid on;
    set(ax, 'GridLineStyle', '--', 'GridAlpha', 0.3);
    set(ax, 'Layer', 'top');  % Grid on top

    % Add legend only to the last subplot
    if m == n_models
        % Reverse order to match visual (Fair on top, Unfair on bottom)
        legend({'Unfair', 'Fair'}, 'Location', 'northeast', ...
               'FontWeight', 'bold', 'Box', 'on');
    end
end

% Adjust layout
set(fig, 'Color', 'w');

% Save combined figure with high resolution
if config.savePNG
    print(gcf, fullfile(resultsDir, 'fm26_individual_fairness_combined.png'), '-dpng', '-r300');
end
if config.savePDF
    print(gcf, fullfile(resultsDir, 'fm26_individual_fairness_combined.pdf'), '-dpdf', '-bestfit');
end
disp("Saved: fm26_individual_fairness_combined.png/pdf");

%% LaTeX Table 2: Individual Fairness Timing (Horizontal Layout)
% Epsilon values as columns, models as rows
disp(" ");
disp("======= INDIVIDUAL FAIRNESS TIMING LATEX TABLE ==========");

latex_timing_filename = fullfile(resultsDir, 'fm26_timing_table.tex');
fid = fopen(latex_timing_filename, 'w');

% Get unique epsilon values for individual fairness (epsilon > 0)
individual_timing = timing_data(timing_data.Epsilon > 0, :);
epsilons_unique = unique(individual_timing.Epsilon);
n_eps = length(epsilons_unique);

fprintf(fid, '\\begin{table}[ht]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Individual Fairness Verification Timing (seconds per sample)}\n');
fprintf(fid, '\\label{tab:individual_timing}\n');

% Create column format with spacing: l for model, then c with padding for each epsilon
col_format = 'l';
for i = 1:n_eps
    col_format = [col_format '@{\hskip 8pt}c'];
end
fprintf(fid, '\\begin{tabular}{%s}\n', col_format);
fprintf(fid, '\\toprule\n');

% Two-row header: first row spans epsilon columns with label
fprintf(fid, ' & \\multicolumn{%d}{c}{Perturbation Level ($\\epsilon$)} \\\\\n', n_eps);
fprintf(fid, '\\cmidrule(l){2-%d}\n', n_eps + 1);

% Second header row with just epsilon values
fprintf(fid, 'Model');
for i = 1:n_eps
    fprintf(fid, ' & %.2f', epsilons_unique(i));
end
fprintf(fid, ' \\\\\n');
fprintf(fid, '\\midrule\n');

% Data rows (one per model) with full display names
for m = 1:length(models)
    modelName = models{m};
    % Use display name if available
    if isKey(model_display_names, modelName)
        displayName = model_display_names(modelName);
    else
        displayName = modelName;
    end
    fprintf(fid, '%s', displayName);

    for i = 1:n_eps
        eps_val = epsilons_unique(i);
        % Find timing for this model and epsilon
        idx = strcmp(individual_timing.Model, modelName) & individual_timing.Epsilon == eps_val;
        if any(idx)
            avg_time = individual_timing.AvgTimePerSample(idx);
            fprintf(fid, ' & %.4f', avg_time);
        else
            fprintf(fid, ' & --');
        end
    end
    fprintf(fid, ' \\\\\n');
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);

% Also print to console
type(latex_timing_filename);
disp(" ");
disp("Saved: " + latex_timing_filename);

%% Summary
disp(" ");
disp("======= FM26 PLOTTING COMPLETE ==========");
disp("Generated outputs in " + resultsDir + ":");
disp("  1. fm26_counterfactual_table.tex (LaTeX table)");
disp("  2. fm26_individual_fairness_combined.png/pdf (Area plot)");
disp("  3. fm26_timing_table.tex (LaTeX table)");
