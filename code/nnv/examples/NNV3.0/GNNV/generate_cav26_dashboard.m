function generate_cav26_dashboard(results, model_data, figuresDir, varargin)
% generate_cav26_dashboard - Create publication-quality dashboard for CAV26
%
% Generates a single dashboard with two sections:
%   - Left: GINE network topology across 3 epsilon values (0.001, 0.005, 0.01)
%   - Right: Epsilon sensitivity chart comparing all models (GCN, GINE, GINE+Edge)
%
% Usage:
%   generate_cav26_dashboard(results, model_data, figuresDir)
%   generate_cav26_dashboard(results, model_data, figuresDir, 'layout', 'manual')
%
% Inputs:
%   results    - Results struct from run_gnn_experiments
%   model_data - Model data struct (from gine_ieee24.mat)
%   figuresDir - Output directory for figures
%
% Options:
%   'layout'  - 'force' (default) or 'manual' for node positioning
%
% Output:
%   figures/dashboard.pdf, dashboard.png
%
% Author: Anne Tumlin
% Date: 01/21/2026

%% Parse options
p = inputParser;
addParameter(p, 'layout', 'force', @(x) ismember(x, {'force', 'manual'}));
parse(p, varargin{:});
opts = p.Results;

%% Color scheme (colorblind-friendly)
colors = struct();
colors.verified = [0.2, 0.65, 0.3];   % Green
colors.unknown  = [0.95, 0.6, 0.1];   % Orange
colors.violated = [0.8, 0.15, 0.15];  % Red
colors.na       = [0.5, 0.5, 0.5];    % Gray

%% Get graph structure from model data
if isfield(model_data, 'src') && isfield(model_data, 'dst')
    src = double(model_data.src);
    dst = double(model_data.dst);
    adj_list = [src, dst];
    numNodes = max([src; dst]);
elseif isfield(model_data, 'ANorm_g')
    A = double(model_data.ANorm_g);
    [src, dst] = find(A);
    adj_list = [src, dst];
    numNodes = size(A, 1);
else
    error('Cannot determine graph structure from model data');
end

%% Get layout coordinates
[x, y] = get_ieee24_layout(opts.layout, adj_list, numNodes);

%% Configuration
models = results.config.models;
epsilons = results.config.epsilons;
num_scenarios = results.config.num_scenarios;

% Total nodes for percentage calculation
total_nodes = results.data{1,1,1}.verified + results.data{1,1,1}.unknown + results.data{1,1,1}.violated;

%% Create figure with manual positioning for clean layout
fig = figure('Position', [100, 100, 1000, 400], 'Visible', 'off');
set(fig, 'DefaultAxesFontName', 'Times New Roman');
set(fig, 'DefaultAxesFontSize', 10);
set(fig, 'DefaultTextFontName', 'Times New Roman');
set(fig, 'Color', 'white');

%% Define layout positions (normalized coordinates [left, bottom, width, height])
% Topology section: left 60% of figure, with margins
topo_left = 0.03;
topo_bottom = 0.18;
topo_width = 0.55;
topo_height = 0.68;
topo_panel_width = topo_width / 3;

% Sensitivity chart: right side - aligned with topology box
sens_left = 0.65;
sens_bottom = 0.10;
sens_width = 0.32;
sens_height = 0.74;

%% Panels 1-3: GINE Network Topology across epsilon values
model_idx = 2;  % GINE model (index 2 in {'GCN', 'GINE', 'GINE+Edge'})
scenario_idx = 1;

% Store handles for shared legend
legend_handles = [];
topo_axes = gobjects(1, 3);

for e = 1:length(epsilons)
    % Create axis for this topology panel
    ax_left = topo_left + (e-1) * topo_panel_width;
    topo_axes(e) = axes('Position', [ax_left, topo_bottom, topo_panel_width, topo_height]);
    hold on;

    % Get verification results for this epsilon value
    exp_data = results.data{model_idx, e, scenario_idx};
    verif_results = exp_data.verif_per_node;

    % Draw edges
    for i = 1:size(adj_list, 1)
        s = adj_list(i, 1);
        t_node = adj_list(i, 2);
        if s <= length(x) && t_node <= length(x)
            plot([x(s), x(t_node)], [y(s), y(t_node)], '-', ...
                'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.8);
        end
    end

    % Draw nodes by status (for legend ordering)
    h_verified = [];
    h_unknown = [];
    h_na = [];

    for i = 1:length(x)
        if i <= length(verif_results)
            status = verif_results(i);
        else
            status = -1;
        end

        switch status
            case 1
                c = colors.verified;
                h = scatter(x(i), y(i), 70, c, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                if isempty(h_verified), h_verified = h; end
            case 2
                c = colors.unknown;
                h = scatter(x(i), y(i), 70, c, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                if isempty(h_unknown), h_unknown = h; end
            case 0
                c = colors.violated;
                scatter(x(i), y(i), 70, c, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            otherwise
                c = colors.na;
                h = scatter(x(i), y(i), 70, c, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                if isempty(h_na), h_na = h; end
        end
    end

    % Store handles from first panel for shared legend
    if e == 1
        legend_handles = [h_verified, h_unknown, h_na];
    end

    % Add node labels
    for i = 1:length(x)
        text(x(i), y(i) - 0.7, sprintf('%d', i), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'FontWeight', 'bold', ...
            'FontName', 'Times New Roman');
    end

    hold off;
    axis equal; axis off;
end

%% Add box around the topology section
box_bottom = topo_bottom - 0.08;
box_height = topo_height + 0.06;
annotation(fig, 'rectangle', [topo_left - 0.01, box_bottom, topo_width + 0.02, box_height], ...
    'LineWidth', 1.5, 'Color', [0.3, 0.3, 0.3]);

%% Add x-axis line with tick marks at bottom of topology section
axis_y = box_bottom + 0.12;
tick_height = 0.015;

% Draw horizontal line
annotation(fig, 'line', [topo_left, topo_left + topo_width], [axis_y, axis_y], ...
    'LineWidth', 1.2, 'Color', 'k');

% Add tick marks and labels for each epsilon
for e = 1:length(epsilons)
    tick_x = topo_left + (e - 0.5) * topo_panel_width;  % Center of each panel
    % Tick mark
    annotation(fig, 'line', [tick_x, tick_x], [axis_y, axis_y - tick_height], ...
        'LineWidth', 1.2, 'Color', 'k');
    % Tick label (epsilon value)
    annotation(fig, 'textbox', [tick_x - 0.04, axis_y - 0.045, 0.08, 0.03], ...
        'String', sprintf('%.3f', epsilons(e)), 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'FontSize', 11, ...
        'FontName', 'Times New Roman');
end

% Add x-axis label below tick labels
annotation(fig, 'textbox', [topo_left, axis_y - 0.095, topo_width, 0.03], ...
    'String', 'Perturbation \epsilon', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');

%% Add titles at same level for both sections
title_y = 0.92;  % Same Y position for both titles

% Topology title with GINE marker (orange line + square to match sensitivity chart)
gine_color = [0.9, 0.6, 0.0];  % Orange - same as GINE in sensitivity chart
title_center_x = topo_left + topo_width/2;  % Center of topology section

% Title text with opening parenthesis - centered
annotation(fig, 'textbox', [topo_left - 0.02, title_y - 0.05, topo_width * 0.52, 0.06], ...
    'String', 'IEEE 24-Bus GINE (', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'right', 'FontSize', 13, 'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');

% Add orange line segment with square marker - centered in title
marker_x = title_center_x - 0.002;
marker_y = title_y - 0.02;
line_half_width = 0.012;
annotation(fig, 'line', [marker_x - line_half_width, marker_x + line_half_width], ...
    [marker_y, marker_y], 'LineWidth', 1.5, 'Color', gine_color);
% Square marker in center
annotation(fig, 'rectangle', [marker_x - 0.004, marker_y - 0.012, 0.008, 0.024], ...
    'FaceColor', gine_color, 'Color', gine_color, 'LineWidth', 1);

% Closing parenthesis and "Verification"
annotation(fig, 'textbox', [marker_x + line_half_width - 0.005, title_y - 0.05, 0.15, 0.06], ...
    'String', ') Verification', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left', 'FontSize', 13, 'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');

%% Add shared legend for node verification status (below topologies)
% Create a dummy axes for the legend
ax_legend = axes('Position', [topo_left, 0.03, topo_width, 0.05], 'Visible', 'off');
hold(ax_legend, 'on');
h1 = scatter(ax_legend, nan, nan, 60, colors.verified, 'filled', 'MarkerEdgeColor', 'k');
h2 = scatter(ax_legend, nan, nan, 60, colors.unknown, 'filled', 'MarkerEdgeColor', 'k');
h3 = scatter(ax_legend, nan, nan, 60, colors.na, 'filled', 'MarkerEdgeColor', 'k');
lgd = legend([h1, h2, h3], {'Verified Safe', 'Unknown', 'Non-Voltage Bus'}, ...
    'Orientation', 'horizontal', 'FontSize', 11, 'FontName', 'Times New Roman', ...
    'Location', 'north');
lgd.Box = 'off';

%% Sensitivity Chart (right side)
ax_sens = axes('Position', [sens_left, sens_bottom, sens_width, sens_height]);
hold on;

% Calculate total percentage: sum verified across all scenarios / (total_nodes * num_scenarios)
total_possible = total_nodes * num_scenarios;  % 13 * 10 = 130
pct_total = zeros(length(models), length(epsilons));

for m = 1:length(models)
    for e = 1:length(epsilons)
        total_verified = 0;
        for s = 1:num_scenarios
            total_verified = total_verified + results.data{m, e, s}.verified;
        end
        pct_total(m, e) = (total_verified / total_possible) * 100;
    end
end

% Plot without error bars (clean lines)
markers = {'o-', 's-', 'd-'};
line_colors = [0.0, 0.45, 0.7;   % Blue (GCN)
               0.9, 0.6, 0.0;    % Orange (GINE)
               0.0, 0.6, 0.5];   % Teal (GINE+Edge)

for m = 1:length(models)
    plot(epsilons, pct_total(m, :), markers{m}, ...
        'Color', line_colors(m, :), 'LineWidth', 2.0, 'MarkerSize', 8, ...
        'MarkerFaceColor', line_colors(m, :));
end
hold off;

xlabel('Perturbation \epsilon', 'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Verified (%)', 'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend(models, 'Location', 'southwest', 'FontSize', 9, 'FontName', 'Times New Roman');

% Add title at same level as topology title (using annotation)
annotation(fig, 'textbox', [sens_left, title_y - 0.05, sens_width, 0.06], ...
    'String', 'Verification vs. Perturbation', 'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold', ...
    'FontName', 'Times New Roman');
xlim([0, max(epsilons) * 1.15]);
ylim([0, 100]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);

%% Save figure
save_figure(fig, fullfile(figuresDir, 'dashboard'));
close(fig);

fprintf('Dashboard saved to: %s\n', figuresDir);
end


%% =========================================================================
%  FIGURE SAVE HELPER
%  =========================================================================

function save_figure(fig, filepath_base)
% Save figure to PDF and PNG

    % Save PNG at moderate resolution for easier viewing
    print(fig, [filepath_base, '.png'], '-dpng', '-r150');

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
%  LAYOUT FUNCTIONS
%  =========================================================================

function [x, y] = get_ieee24_layout(method, adj_list, numNodes)
% Get x,y coordinates for IEEE 24-bus system

    switch method
        case 'force'
            [x, y] = force_directed_layout(adj_list, numNodes);
        case 'manual'
            [x, y] = manual_ieee24_layout();
        otherwise
            error('Unknown layout method: %s', method);
    end
end


function [x, y] = force_directed_layout(adj_list, numNodes)
% Compute force-directed layout using MATLAB's graph functions

    G = graph(adj_list(:,1), adj_list(:,2));
    fig_temp = figure('Visible', 'off');
    h = plot(G, 'Layout', 'force', 'Iterations', 200);
    x = h.XData';
    y = h.YData';
    close(fig_temp);

    % Normalize to [0, 12] range
    x = (x - min(x)) / (max(x) - min(x)) * 12;
    y = (y - min(y)) / (max(y) - min(y)) * 12;
end


function [x, y] = manual_ieee24_layout()
% Hand-tuned coordinates for IEEE 24-bus RTS

    x = zeros(24, 1);
    y = zeros(24, 1);

    % Area 1: Buses 1-10
    x(1) = 1; y(1) = 1;
    x(2) = 2; y(2) = 1;
    x(3) = 3; y(3) = 1;
    x(4) = 2; y(4) = 2;
    x(5) = 3; y(5) = 2;
    x(6) = 1; y(6) = 3;
    x(7) = 2; y(7) = 3;
    x(8) = 3; y(8) = 3;
    x(9) = 1.5; y(9) = 4;
    x(10) = 2.5; y(10) = 4;

    % Area 2: Buses 11-24
    x(11) = 5; y(11) = 3;
    x(12) = 5; y(12) = 4;
    x(13) = 6; y(13) = 1;
    x(14) = 7; y(14) = 1;
    x(15) = 8; y(15) = 1;
    x(16) = 6; y(16) = 2;
    x(17) = 7; y(17) = 2;
    x(18) = 8; y(18) = 2;
    x(19) = 6; y(19) = 3;
    x(20) = 7; y(20) = 3;
    x(21) = 8; y(21) = 3;
    x(22) = 6; y(22) = 4;
    x(23) = 7; y(23) = 4;
    x(24) = 8; y(24) = 4;

    % Scale
    x = x * 1.1;
    y = y * 2;
end
