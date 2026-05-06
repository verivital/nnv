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
    % TLLverify dropped from active comparison; rows still in
    % acas_rl_tll/legacy/results_tllverify.{mat,csv}.
    benches = {'acas_p3','acas_p4','rl','oval21','collins_rul'};

    % Per-benchmark headline algorithm grids (matches the driver's
    % algorithms_for(tool, benchmark) defaults: full grid for ACAS,
    % approx-star only for RL/oval21/collins_rul). Anything else in
    % the .mat is treated as "supplement" data and excluded from the
    % headline table.
    headline_algs = containers.Map();
    headline_algs('acas_p3') = {'approx-star','relax-star-range-25', ...
        'relax-star-range-50','relax-star-range-75','relax-star-range-100', ...
        'exact-star','estimate-bounds','alpha-beta-crown'};
    headline_algs('acas_p4')      = headline_algs('acas_p3');
    headline_algs('rl')           = {'approx-star','relax-star-range-50', ...
        'exact-star','estimate-bounds'};
    headline_algs('oval21')       = {'approx-star', ...
        'relax-star-area-25','relax-star-area-50', ...
        'relax-star-area-75','relax-star-area-100','estimate-bounds'};
    headline_algs('collins_rul')  = {'approx-star', ...
        'relax-star-area-25','relax-star-area-50', ...
        'relax-star-area-75','relax-star-area-100', ...
        'exact-star','estimate-bounds'};

    % Column order: counts, then per-row timeout (read from the .mat row's
    % `timeout` column to expose mixed regimes — 300s for ACAS/RL global cap
    % vs per-instance VNNCOMP timeouts on OVAL21/Collins RUL), then mean
    % wall-time on solved instances, then PAR-2 (VNN-COMP scoring metric:
    % unsolved instances counted at 2*timeout). PAR-2 is what NN-verification
    % reviewers expect; without it any quoted "mean time" is misleading
    % because timeouts get filtered out of the average.
    header = {"Benchmark","Tool","Algorithm","V","X","?","T/O","Timeout (s)","Mean t (s)","PAR-2 (s)"};
    rows = {};
    suppRows = {};
    txtLines = strings(0,1);
    for i = 1:numel(benches)
        bench   = benches{i};
        matFile = fullfile(opts.resultsDir, sprintf("results_%s.mat", bench));
        R = u.load_results(matFile);
        if isempty(R)
            txtLines(end+1,1) = sprintf("%-12s no results yet", bench); %#ok<AGROW>
            continue;
        end
        algSet = headline_algs(bench);
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
            % Per-instance timeout: most rows in a (tool, alg, bench) group
            % share a timeout. Surface the mode; if it varies (per-instance
            % VNNCOMP caps), append "*" to flag it.
            timeouts = sel.timeout;
            if numel(unique(timeouts)) <= 1
                tStr = u.format_time(timeouts(1), Inf);
            else
                tStr = sprintf("%s*", u.format_time(median(timeouts), Inf));
            end
            % PAR-2 score: unsolved instances counted at 2*timeout.
            % Use median timeout when per-instance timeouts vary.
            tEff = median(timeouts(~isnan(timeouts)));
            if isnan(tEff), tEff = 0; end
            par2 = u.par2(sel.time, sel.status, tEff);
            rowCells = { bench, tool, alg, ...
                sprintf("%d", v), sprintf("%d", vi), sprintf("%d", un), sprintf("%d", to), ...
                tStr, u.format_time(mt), u.format_time(par2) };
            isHeadline = any(string(algSet) == alg);
            if isHeadline
                rows{end+1} = rowCells; %#ok<AGROW>
                txtLines(end+1,1) = sprintf("%-12s %-14s %-22s V=%3d X=%3d ?=%3d T/O=%3d  T=%s  mean=%s  PAR-2=%s", ...
                    bench, tool, alg, v, vi, un, to, tStr, u.format_time(mt), u.format_time(par2)); %#ok<AGROW>
            else
                suppRows{end+1} = rowCells; %#ok<AGROW>
            end
        end
    end
    if ~isempty(suppRows)
        txtLines(end+1,1) = "";
        txtLines(end+1,1) = "--- supplement (out-of-headline-grid rows; preserved in .mat) ---";
        for k = 1:numel(suppRows)
            r = suppRows{k};
            txtLines(end+1,1) = sprintf("%-12s %-14s %-22s V=%3s X=%3s ?=%3s T/O=%3s  T=%s  mean=%s  PAR-2=%s", ...
                r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, r{7}, r{8}, r{9}, r{10}); %#ok<AGROW>
        end
    end
    txtLines(end+1,1) = "";
    txtLines(end+1,1) = "Timeout column: '*' suffix = per-instance timeouts vary (VNNCOMP-style); value shown is the median.";

    u.emit_latex_table( fullfile(opts.outDir,'table_A.tex'), ...
        header, rows, ...
        "ToolComparison FC-net half: NNV vs MathWorks AIVL on ACAS Xu and RL controllers.", ...
        "tab:toolcomparison-acas-rl");
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
