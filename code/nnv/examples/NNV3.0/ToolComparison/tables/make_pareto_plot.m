function make_pareto_plot(varargin)
%MAKE_PARETO_PLOT Per-benchmark Pareto plot of verification rate vs PAR-2.
%
%   For each benchmark in acas_rl_tll/results/results_<benchmark>.mat,
%   plots one point per (tool, algorithm) at:
%     x = PAR-2 score (s)              -- "cost"  (lower is better)
%     y = verification rate (V/total)  -- "value" (higher is better)
%
%   The Pareto frontier (lower-right is best: low PAR-2 + high V-rate) is
%   highlighted; dominated points are dimmed. Tool is encoded by marker
%   color (nnv=blue, aivl=orange) and algorithm by marker shape.
%
%   Outputs:
%     tables/out/pareto_<benchmark>.png
%     tables/out/pareto_<benchmark>.pdf
%
%   Additive companion to the headline tables — does NOT replace them per
%   PI directive. Intended as the at-a-glance figure for paper / talks.
%
%   Caveat: a "Pareto" view assumes (V-rate, PAR-2) are the right two axes.
%   This loses violation information (X count) — for benchmarks where a
%   tool finds counterexamples (e.g., RL controllers, OVAL21), high V-rate
%   alone undersells AIVL's contribution. Read alongside Table A.

    p = inputParser;
    addParameter(p, 'resultsDir', ...
        fullfile(toolcomparison_root(), 'acas_rl_tll', 'results'));
    addParameter(p, 'outDir', fullfile(toolcomparison_root(), 'tables', 'out'));
    addParameter(p, 'benchmarks', {'acas_p3','acas_p4','rl','oval21','collins_rul'});
    addParameter(p, 'show', false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;
    if ~isfolder(opts.outDir), mkdir(opts.outDir); end

    u = tool_utils();

    for i = 1:numel(opts.benchmarks)
        bench   = opts.benchmarks{i};
        matFile = fullfile(opts.resultsDir, sprintf("results_%s.mat", bench));
        R = u.load_results(matFile);
        if isempty(R), continue; end
        plot_one(bench, R, u, opts);
    end
end

% =========================================================================

function plot_one(bench, R, u, opts)
    keys = unique(R(:, {'tool','algorithm'}));
    n = height(keys);
    if n == 0, return; end

    pts = struct('tool', {}, 'algorithm', {}, 'vrate', {}, 'par2', {});
    for k = 1:n
        tool = keys.tool(k);
        alg  = keys.algorithm(k);
        sel  = R(R.tool==tool & R.algorithm==alg & R.benchmark==string(bench), :);
        if height(sel) == 0, continue; end
        vrate = sum(sel.status == "verified") / height(sel);
        tEff  = median(sel.timeout(~isnan(sel.timeout)));
        if isnan(tEff), tEff = 0; end
        par2  = u.par2(sel.time, sel.status, tEff);
        pts(end+1) = struct('tool', tool, 'algorithm', alg, ...
                            'vrate', vrate, 'par2', par2); %#ok<AGROW>
    end
    if isempty(pts), return; end

    [is_pareto, dominated] = compute_pareto(pts);

    fig = figure('Visible', tern(opts.show, 'on', 'off'), ...
                 'Position', [100 100 700 500]);
    hold on; grid on;
    set(gca, 'XScale', 'log');

    % Plot dominated points dimmed.
    for k = dominated
        plot(pts(k).par2, pts(k).vrate, marker_for(pts(k)), ...
             'MarkerSize', 8, 'MarkerFaceColor', color_for(pts(k).tool, 0.3), ...
             'MarkerEdgeColor', [0.6 0.6 0.6], 'LineWidth', 0.5);
    end
    % Plot Pareto-frontier points solid.
    for k = is_pareto
        plot(pts(k).par2, pts(k).vrate, marker_for(pts(k)), ...
             'MarkerSize', 12, 'MarkerFaceColor', color_for(pts(k).tool, 1.0), ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
    end
    % Connect Pareto frontier in increasing PAR-2 order.
    if numel(is_pareto) > 1
        front = pts(is_pareto);
        [~, ord] = sort([front.par2]);
        front = front(ord);
        plot([front.par2], [front.vrate], '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 0.7);
    end

    % Per-point labels.
    for k = 1:numel(pts)
        text(pts(k).par2, pts(k).vrate + 0.025, ...
             sprintf("%s/%s", pts(k).tool, short_alg(pts(k).algorithm)), ...
             'FontSize', 7, 'HorizontalAlignment', 'center', ...
             'Interpreter', 'none');
    end

    xlabel("PAR-2 score (s)  \rightarrow  lower is better");
    ylabel("Verification rate (V / total)  \rightarrow  higher is better");
    title(sprintf("ToolComparison Pareto: %s", bench), 'Interpreter', 'none');
    ylim([-0.05, 1.05]);

    base = fullfile(opts.outDir, sprintf("pareto_%s", bench));
    exportgraphics(fig, base + ".png", 'Resolution', 200);
    exportgraphics(fig, base + ".pdf", 'ContentType', 'vector');
    if ~opts.show, close(fig); end
    fprintf("Wrote %s.{png,pdf}\n", base);
end

function [is_pareto, dominated] = compute_pareto(pts)
% A point is on the Pareto frontier if no other point has BOTH lower or
% equal par2 AND higher or equal vrate (with at least one strict).
    n = numel(pts);
    is_pareto = false(1, n);
    for k = 1:n
        dominated_by_other = false;
        for j = 1:n
            if j == k, continue; end
            if pts(j).par2 <= pts(k).par2 && pts(j).vrate >= pts(k).vrate && ...
               (pts(j).par2 < pts(k).par2 || pts(j).vrate > pts(k).vrate)
                dominated_by_other = true; break;
            end
        end
        is_pareto(k) = ~dominated_by_other;
    end
    dominated = find(~is_pareto);
    is_pareto = find(is_pareto);
end

function c = color_for(tool, alpha)
    switch string(tool)
        case "nnv",  base = [0.0   0.4470 0.7410];     % blue
        case "aivl", base = [0.85  0.3250 0.0980];     % orange-red
        otherwise,   base = [0.5   0.5    0.5];        % gray
    end
    % Apply alpha as desaturation.
    c = base * alpha + (1 - alpha) * [1 1 1];
end

function m = marker_for(pt)
% Map algorithm to a stable marker. Falls back to 'o'.
    switch string(pt.algorithm)
        case "approx-star",         m = 'o';   % circle
        case "exact-star",          m = 's';   % square
        case "estimate-bounds",     m = 'd';   % diamond
        case "deep-poly",           m = 'p';   % pentagon
        case "alpha-beta-crown",    m = 'h';   % hexagon
        case {"relax-star-range-25","relax-star-area-25"},   m = '^';  % up triangle
        case {"relax-star-range-50","relax-star-area-50"},   m = 'v';  % down triangle
        case {"relax-star-range-75","relax-star-area-75"},   m = '<';  % left triangle
        case {"relax-star-range-100","relax-star-area-100"}, m = '>';  % right triangle
        otherwise, m = 'o';
    end
end

function s = short_alg(alg)
% Compact algorithm label for in-plot text.
    map = {
        "approx-star",            "approx";
        "exact-star",             "exact";
        "relax-star-range-25",    "rlx-r-25";
        "relax-star-range-50",    "rlx-r-50";
        "relax-star-range-75",    "rlx-r-75";
        "relax-star-range-100",   "rlx-r-100";
        "relax-star-area-25",     "rlx-a-25";
        "relax-star-area-50",     "rlx-a-50";
        "relax-star-area-75",     "rlx-a-75";
        "relax-star-area-100",    "rlx-a-100";
        "estimate-bounds",        "est-bnds";
        "deep-poly",              "deep-poly";
        "alpha-beta-crown",       "ab-CROWN";
    };
    s = string(alg);
    for k = 1:size(map, 1)
        if s == map{k, 1}, s = map{k, 2}; return; end
    end
end

function y = tern(c, a, b), if c, y = a; else, y = b; end, end

function r = toolcomparison_root()
    r = fileparts(fileparts(mfilename('fullpath')));
end
