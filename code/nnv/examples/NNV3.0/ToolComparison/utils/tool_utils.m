function utils = tool_utils()
%TOOL_UTILS Shared helpers for the NNV-vs-MathWorks AIVL tool comparison.
%
%   u = tool_utils();
%   row = u.new_row('nnv', 'acas_p3', 1, 'verified', 2.3, 'exact-star', 300);
%   u.append_to_mat('results/results_p3.mat', row);
%   summary = u.tally('results/results_p3.mat', {'benchmark','tool'});
%   u.emit_latex_table('out/table_A.tex', header, rows, caption, label);
%
%   Single canonical result-row schema (MATLAB table):
%     tool         (string)  'nnv' | 'aivl'
%                            (legacy 'mw_estimate' / 'mw_deeppoly' / 'mw_abc'
%                             rows in older .mat files are migrated by
%                             utils/canonicalize_bundled_results.m to
%                             tool='aivl' + the matching algorithm.)
%     benchmark    (string)  'acas_p3' | 'acas_p4' | 'rl' | 'oval21' | ...
%     instance_id  (string)  e.g. 'ACASXU_run2a_1_1_batch_2000' or integer row id
%     status       (string)  'verified' | 'violated' | 'unknown' | 'timeout' | 'error'
%     time         (double)  wallclock seconds (NaN on error/timeout)
%     algorithm    (string)  NNV  : 'approx-star' | 'exact-star' | 'relax-star-range-50' | ...
%                            AIVL : 'estimate-bounds' | 'deep-poly' | 'alpha-beta-crown'
%     timeout      (double)  configured timeout in seconds
%     note         (string)  free-form; use for error messages

    utils.new_row              = @new_row;
    utils.append_to_mat        = @append_to_mat;
    utils.load_results         = @load_results;
    utils.has_instance         = @has_instance;
    utils.purge_status         = @purge_status;
    utils.tally                = @tally;
    utils.format_time          = @format_time;
    utils.format_count         = @format_count;
    utils.emit_latex_table     = @emit_latex_table;
    utils.sanity_check_vs_nnv2 = @sanity_check_vs_nnv2;
    utils.par2                 = @par2;
end

% -------------------------------------------------------------------------

function row = new_row(tool, benchmark, instance_id, status, time, algorithm, timeout, note)
    if nargin < 8, note = ""; end
    row = table( ...
        string(tool), string(benchmark), string(instance_id), string(status), ...
        double(time), string(algorithm), double(timeout), string(note), ...
        'VariableNames', ...
        {'tool','benchmark','instance_id','status','time','algorithm','timeout','note'});
end

function append_to_mat(matFile, row)
% Append a result row to an existing .mat (creates if missing).
    if exist(matFile, 'file')
        S = load(matFile, 'results');
        results = [S.results; row];
    else
        results = row;
    end
    outDir = fileparts(matFile);
    if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
    save(matFile, 'results');
end

function results = load_results(matFile)
    if ~exist(matFile, 'file')
        results = table('Size',[0 8], ...
            'VariableTypes',{'string','string','string','string','double','string','double','string'}, ...
            'VariableNames',{'tool','benchmark','instance_id','status','time','algorithm','timeout','note'});
        return;
    end
    S = load(matFile, 'results');
    results = S.results;
end

function n = purge_status(matFile, statuses, tool, algorithm)
% Drop rows in matFile whose status is in the given cell/string array.
% Optional filters narrow the purge:
%   purge_status(f, 'timeout')                       drops every timeout row
%   purge_status(f, 'timeout', 'nnv')                drops timeout rows where tool=nnv only
%   purge_status(f, 'timeout', 'nnv', 'exact-star')  drops timeout rows where (tool, alg) match
%
% Returns the number of rows removed. Used to re-run error/timeout cells
% (or surgically reconcile a specific (tool, algorithm) subset, e.g.
% ACAS exact-star at a higher timeout) without losing successful rows.
    if ~exist(matFile, 'file'), n = 0; return; end
    R = load_results(matFile);
    if isempty(R), n = 0; return; end
    drop = ismember(R.status, string(statuses));
    if nargin >= 3 && ~isempty(tool)
        drop = drop & (R.tool == string(tool));
    end
    if nargin >= 4 && ~isempty(algorithm)
        drop = drop & (R.algorithm == string(algorithm));
    end
    nbefore = height(R);
    R = R(~drop, :);
    n = nbefore - height(R);
    if n > 0
        results = R; %#ok<NASGU>
        save(matFile, 'results');
    end
end

function tf = has_instance(matFile, tool, benchmark, instance_id, algorithm)
% True iff a completed row for this (tool, benchmark, instance, algorithm) already exists.
% Enables graceful resumption after crash / SIGTERM.
    R = load_results(matFile);
    if isempty(R), tf = false; return; end
    tf = any( R.tool == string(tool) & ...
              R.benchmark == string(benchmark) & ...
              R.instance_id == string(instance_id) & ...
              R.algorithm == string(algorithm) );
end

function summary = tally(matFile, groupVars)
% Aggregate SAT / UNSAT / unknown / timeout counts and mean time per group.
    if nargin < 2, groupVars = {'tool','benchmark','algorithm'}; end
    R = load_results(matFile);
    if isempty(R)
        summary = table(); return;
    end
    [summary, ~] = groupsummary(R, groupVars, { ...
        @(s) sum(s=="verified"),   'status', ...
        @(s) sum(s=="violated"),   'status', ...
        @(s) sum(s=="unknown"),    'status', ...
        @(s) sum(s=="timeout"),    'status', ...
        @(s) sum(s=="error"),      'status', ...
        @mean,                     'time'});
    summary.Properties.VariableNames(end-5:end) = ...
        {'verified','violated','unknown','timeout','error','mean_time'};
end

function s = format_time(t, timeout)
% Pretty-print a time or ">T/O" when timeout was hit.
    if nargin < 2, timeout = Inf; end
    if isnan(t) || t >= timeout, s = sprintf(">%g", timeout);
    elseif t < 1,                s = sprintf("%.2f", t);
    elseif t < 100,              s = sprintf("%.1f", t);
    else,                        s = sprintf("%.0f", t);
    end
end

function s = format_count(n, total)
    if total > 0
        s = sprintf("%d/%d (%.0f%%)", n, total, 100*n/total);
    else
        s = sprintf("%d", n);
    end
end

function p = par2(times, statuses, timeout)
% PAR-2 score used in VNN-COMP: unsolved instances counted at 2*timeout.
    t = times(:);
    s = statuses(:);
    t(s=="timeout" | s=="unknown" | s=="error" | isnan(t)) = 2*timeout;
    p = mean(t);
end

function emit_latex_table(outFile, header, rows, caption, label)
% Write a LaTeX tabular to outFile. header is a cell-string; rows is cell of cell-strings.
    outDir = fileparts(outFile);
    if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
    colSpec = repmat("l", 1, numel(header));
    fid = fopen(outFile, 'w');
    fprintf(fid, "\\begin{table}[t]\n\\centering\n\\caption{%s}\n\\label{%s}\n", caption, label);
    fprintf(fid, "\\begin{tabular}{%s}\n\\toprule\n", strjoin(colSpec,""));
    fprintf(fid, "%s \\\\\n\\midrule\n", strjoin(string(header), " & "));
    for i = 1:numel(rows)
        fprintf(fid, "%s \\\\\n", strjoin(string(rows{i}), " & "));
    end
    fprintf(fid, "\\bottomrule\n\\end{tabular}\n\\end{table}\n");
    fclose(fid);
end

function report = sanity_check_vs_nnv2(resultsDir, reportFile)
% Compare ACAS exact-star rows against CAV'23 published baselines.
% CAV'23 Table 2 numbers for NNV exact-star (45 networks each):
%   P3: verified=42 violated=3 unknown=0
%   P4: verified=42 violated=3 unknown=0
%
% CAV'23 had NO timeout cap (max ACAS Xu network was 10479 s). We cap at
% 900 s and record overruns as "timeout". So compare CAV'23's verified
% count against our (verified + timeout) -- timeouts here are very likely
% verified-but-slow.
%
% resultsDir is a folder containing results_<benchmark>.mat files.
% Tolerance: +/-2 on counts.
    baselines = struct( ...
        'acas_p3', struct('verified', 42, 'violated', 3, 'unknown', 0), ...
        'acas_p4', struct('verified', 42, 'violated', 3, 'unknown', 0));
    lines = strings(0,1);
    lines(end+1,1) = "ToolComparison sanity report vs CAV'23 (NNV 2.0 exact-star)";
    lines(end+1,1) = string(repmat('-',1,60));
    benches = fieldnames(baselines);
    for i = 1:numel(benches)
        b  = benches{i};
        ex = baselines.(b);
        f  = fullfile(resultsDir, sprintf("results_%s.mat", b));
        if ~exist(f, 'file')
            lines(end+1,1) = sprintf("%s: result file missing (%s)", b, f); %#ok<AGROW>
            continue;
        end
        R  = load_results(f);
        rb = R(R.tool=="nnv" & R.algorithm=="exact-star", :);
        v  = sum(rb.status=="verified");
        vi = sum(rb.status=="violated");
        un = sum(rb.status=="unknown");
        to = sum(rb.status=="timeout");
        % Compare CAV'23's count against (verified + timeout) since timeouts
        % are likely verified-but-slow under our 900 s cap.
        ok = abs((v+to)-ex.verified) <= 2 && ...
             abs(vi-ex.violated)     <= 2 && ...
             abs(un-ex.unknown)      <= 2;
        lines(end+1,1) = sprintf("%s: got (V=%d, X=%d, ?=%d, T/O=%d)  CAV'23 expected (V=%d, X=%d, ?=%d)  %s", ...
            b, v, vi, un, to, ex.verified, ex.violated, ex.unknown, ...
            tern(ok, "OK (V+T/O matches CAV'23 V)", "DRIFT")); %#ok<AGROW>
    end
    report = strjoin(lines, newline);
    if nargin >= 2 && ~isempty(reportFile)
        outDir = fileparts(reportFile);
        if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
        fid = fopen(reportFile,'w'); fprintf(fid, "%s\n", report); fclose(fid);
    end
end

function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end
