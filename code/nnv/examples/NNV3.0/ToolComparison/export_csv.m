function export_csv(varargin)
%EXPORT_CSV Export every result .mat in the ToolComparison tree to a CSV alongside.
%
%   export_csv() walks acas_rl_tll/results/*.mat and mnist_resnet/results/*.mat,
%   writes the canonical result table to a CSV with the same basename, and
%   also produces *_summary.csv -- one row per (tool, algorithm, status)
%   group with counts and mean/median time.
%
%   The CSVs open natively in VS Code (CSV preview), Excel, pandas, etc.
%
%   Per-instance file e.g. results_acas_p3.csv:
%     tool,benchmark,instance_id,status,time,algorithm,timeout,note
%
%   Summary file e.g. results_acas_p3_summary.csv:
%     tool,algorithm,status,count,mean_time_s,median_time_s,max_time_s

    p = inputParser;
    addParameter(p, 'roots', { ...
        fullfile(fileparts(mfilename('fullpath')), 'acas_rl_tll', 'results'), ...
        fullfile(fileparts(mfilename('fullpath')), 'mnist_resnet', 'results')});
    parse(p, varargin{:});
    opts = p.Results;

    u = tool_utils();
    nfiles = 0;
    for i = 1:numel(opts.roots)
        d = opts.roots{i};
        if ~isfolder(d), continue; end
        mats = dir(fullfile(d, '*.mat'));
        for k = 1:numel(mats)
            matFile = fullfile(mats(k).folder, mats(k).name);
            R = u.load_results(matFile);
            if isempty(R), continue; end
            base = fullfile(mats(k).folder, mats(k).name(1:end-4));
            csvFile = base + ".csv";
            writetable(R, csvFile);
            fprintf("wrote %s (%d rows)\n", csvFile, height(R));

            % Summary
            keys = unique(R(:, {'tool','algorithm','status'}));
            nrows = height(keys);
            S = table('Size',[nrows 6], ...
                'VariableTypes', {'string','string','string','double','double','double'}, ...
                'VariableNames', {'tool','algorithm','status','count','mean_time_s','median_time_s'});
            for j = 1:nrows
                sel = R( R.tool      == keys.tool(j) & ...
                         R.algorithm == keys.algorithm(j) & ...
                         R.status    == keys.status(j), :);
                S.tool(j)         = keys.tool(j);
                S.algorithm(j)    = keys.algorithm(j);
                S.status(j)       = keys.status(j);
                S.count(j)        = height(sel);
                t = sel.time(~isnan(sel.time));
                if isempty(t), S.mean_time_s(j) = NaN; S.median_time_s(j) = NaN;
                else,         S.mean_time_s(j) = mean(t); S.median_time_s(j) = median(t);
                end
            end
            sumFile = base + "_summary.csv";
            writetable(S, sumFile);
            fprintf("wrote %s (%d groups)\n", sumFile, height(S));
            nfiles = nfiles + 1;
        end
    end
    fprintf("\n%d result files exported.\n", nfiles);
end
