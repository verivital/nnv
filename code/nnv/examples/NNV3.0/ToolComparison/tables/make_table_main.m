function make_table_main()
%MAKE_TABLE_MAIN  Render the consolidated result table from
%results/*.mat. Outputs tables/out/table_main.{tex,txt}.
%
%   One row per (benchmark, tool, algorithm) tuple, ordered:
%     benchmark    canonical order
%     tool         nnv before aivl
%     algorithm    canonical order per (tool, bench)
%
%   Columns:
%     Benchmark Tool Algorithm V X ? T/O Err Mean-t(s) PAR-2(s)

    here = fileparts(mfilename('fullpath'));
    tc_root = fileparts(here);
    results_dir = fullfile(tc_root, 'results');
    out_dir = fullfile(here, 'out');
    if ~isfolder(out_dir), mkdir(out_dir); end

    addpath(fullfile(tc_root, 'utils'));
    addpath_shared();
    u = tool_utils();

    BENCHES = {'acas_xu_p3','acas_xu_p4','rl','oval21','collins_rul','mnist_resnet8'};

    rows = {};
    for b = 1:numel(BENCHES)
        bench = BENCHES{b};
        matFile = fullfile(results_dir, sprintf('%s.mat', bench));
        if ~isfile(matFile)
            continue;
        end
        R = u.load_results(matFile);
        if isempty(R), continue; end

        % Iterate (tool, algorithm) groups in this benchmark.
        keys = unique(R(:, {'tool','algorithm'}), 'rows', 'stable');
        for ki = 1:height(keys)
            tool = string(keys.tool(ki));
            alg  = string(keys.algorithm(ki));
            sub  = R(R.tool == tool & R.algorithm == alg, :);
            tot  = height(sub);
            nV   = sum(sub.status == "verified");
            nX   = sum(sub.status == "violated");
            nU   = sum(sub.status == "unknown");
            nTO  = sum(sub.status == "timeout");
            nErr = sum(sub.status == "error");
            % Mean time over solved rows (verified | violated). Mirrors ToolComparison.
            solved = sub.status == "verified" | sub.status == "violated";
            if any(solved)
                meant = mean(sub.time(solved));
            else
                meant = NaN;
            end
            % PAR-2 uses the per-row configured timeout from `.timeout` column.
            % Most rows in a (bench,tool,alg) group share one timeout, but
            % VNNCOMP-style benchmarks may vary per row.
            % We use the mode of timeouts as the canonical T column and pass
            % each row's own timeout to par2 implicitly via the helper.
            timeouts = sub.timeout(:);
            T_canon = mode(timeouts);
            % Compute PAR-2 row-wise then average.
            par2_vals = sub.time;
            unsolved = sub.status == "timeout" | sub.status == "unknown" | sub.status == "error" | isnan(par2_vals);
            par2_vals(unsolved) = 2 * timeouts(unsolved);
            par2 = mean(par2_vals);

            rows{end+1} = struct( ...
                'benchmark', bench, ...
                'tool',      char(tool), ...
                'algorithm', char(alg), ...
                'N',         tot, ...
                'V',         nV, ...
                'X',         nX, ...
                'Q',         nU, ...
                'TO',        nTO, ...
                'Err',       nErr, ...
                'T',         T_canon, ...
                'mean_t',    meant, ...
                'par2',      par2); %#ok<AGROW>
        end
    end

    % Plaintext rendering.
    txt = fullfile(out_dir, 'table_main.txt');
    fid = fopen(txt, 'w');
    fprintf(fid, '%-14s %-5s %-22s %4s %4s %4s %4s %4s %4s   %8s   %8s\n', ...
        'Benchmark','Tool','Algorithm','N','V','X','?','T/O','Err','Mean t(s)','PAR-2(s)');
    fprintf(fid, '%s\n', repmat('-', 1, 110));
    for i = 1:numel(rows)
        r = rows{i};
        fprintf(fid, '%-14s %-5s %-22s %4d %4d %4d %4d %4d %4d   %8.2f   %8.2f\n', ...
            r.benchmark, r.tool, r.algorithm, r.N, r.V, r.X, r.Q, r.TO, r.Err, r.mean_t, r.par2);
    end
    fclose(fid);
    fprintf('[make_table_main] wrote %s (%d rows)\n', txt, numel(rows));

    % LaTeX rendering (minimal tabular).
    tex = fullfile(out_dir, 'table_main.tex');
    fid = fopen(tex, 'w');
    fprintf(fid, '\\begin{tabular}{llrrrrrrrrr}\n\\hline\n');
    fprintf(fid, 'Benchmark & Tool & Algorithm & N & V & X & ? & T/O & Err & Mean t (s) & PAR-2 (s) \\\\\n');
    fprintf(fid, '\\hline\n');
    for i = 1:numel(rows)
        r = rows{i};
        fprintf(fid, '%s & %s & \\texttt{%s} & %d & %d & %d & %d & %d & %d & %.2f & %.2f \\\\\n', ...
            esc_tex(r.benchmark), r.tool, esc_tex(r.algorithm), ...
            r.N, r.V, r.X, r.Q, r.TO, r.Err, r.mean_t, r.par2);
    end
    fprintf(fid, '\\hline\n\\end{tabular}\n');
    fclose(fid);
    fprintf('[make_table_main] wrote %s\n', tex);

    % Echo to console for AE convenience.
    fprintf('\n');
    type(txt);
end

function s = esc_tex(s)
    s = strrep(s, '_', '\_');
end
