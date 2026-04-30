function make_acas_rl_tll_table(varargin)
%MAKE_ACAS_RL_TLL_TABLE Emit Table A (ToolComparison FC-net half).
%
%   Reads ToolComparison/acas_rl_tll/results/results_<benchmark>.mat and writes:
%     out/table_A.tex   LaTeX tabular
%     out/table_A.txt   plain-text mirror
%     out/sanity_report.txt   CAV'23 cross-check (NNV exact-star vs published)
%
%   Tabular layout: per (benchmark, tool, algorithm), counts and mean time:
%     V (verified) / X (violated) / ? (unknown) / T/O (timeout) / mean t (s)

    p = inputParser;
    addParameter(p, 'resultsDir', ...
        fullfile(toolcomparison_root(), 'acas_rl_tll', 'results'));
    addParameter(p, 'outDir', fullfile(toolcomparison_root(), 'tables', 'out'));
    parse(p, varargin{:});
    opts = p.Results;
    if ~isfolder(opts.outDir), mkdir(opts.outDir); end

    u = tool_utils();
    benches = {'acas_p3','acas_p4','rl','tllverify'};

    header = {"Benchmark","Tool","Algorithm","V","X","?","T/O","Mean t (s)"};
    rows = {};
    txtLines = strings(0,1);
    for i = 1:numel(benches)
        bench   = benches{i};
        matFile = fullfile(opts.resultsDir, sprintf("results_%s.mat", bench));
        R = u.load_results(matFile);
        if isempty(R)
            txtLines(end+1,1) = sprintf("%-12s no results yet", bench); %#ok<AGROW>
            continue;
        end
        keys = unique(R(:, {'tool','algorithm'}));
        for k = 1:height(keys)
            tool = keys.tool(k);
            alg  = keys.algorithm(k);
            sel  = R(R.tool==tool & R.algorithm==alg & R.benchmark==string(bench), :);
            v   = sum(sel.status == "verified");
            vi  = sum(sel.status == "violated");
            un  = sum(sel.status == "unknown");
            to  = sum(sel.status == "timeout");
            mt  = mean(sel.time(sel.status ~= "timeout" & sel.status ~= "error" & ~isnan(sel.time)));
            rows{end+1} = { bench, tool, alg, ...
                sprintf("%d", v), sprintf("%d", vi), sprintf("%d", un), sprintf("%d", to), ...
                u.format_time(mt) }; %#ok<AGROW>
            txtLines(end+1,1) = sprintf("%-12s %-14s %-22s V=%3d X=%3d ?=%3d T/O=%3d mean=%s", ...
                bench, tool, alg, v, vi, un, to, u.format_time(mt)); %#ok<AGROW>
        end
    end

    u.emit_latex_table( fullfile(opts.outDir,'table_A.tex'), ...
        header, rows, ...
        "ToolComparison FC-net half: NNV vs MathWorks AIVL on ACAS Xu, RL, TLL.", ...
        "tab:toolcomparison-acas-rl-tll");
    fid = fopen(fullfile(opts.outDir,'table_A.txt'),'w');
    fprintf(fid, "%s\n", strjoin(txtLines, newline));
    fclose(fid);
    fprintf("Wrote %s and %s\n", ...
        fullfile(opts.outDir,'table_A.tex'), ...
        fullfile(opts.outDir,'table_A.txt'));

    % CAV'23 cross-check (only meaningful when exact-star ACAS rows exist).
    sanityFile = fullfile(opts.outDir, 'sanity_report.txt');
    u.sanity_check_vs_nnv2(opts.resultsDir, sanityFile);
    fprintf("Wrote %s\n", sanityFile);
end

function r = toolcomparison_root()
    % This file: ToolComparison/tables/make_acas_rl_tll_table.m
    %   -> peel filename -> .../tables/  -> peel  -> .../ToolComparison/
    r = fileparts(fileparts(mfilename('fullpath')));
end
