function make_mnist_resnet_table(varargin)
%MAKE_MNIST_RESNET_TABLE Emit Table C (ToolComparison ResNet head-to-head).
%
%   Reads ToolComparison/mnist_resnet/results/expC_<model>.mat and writes:
%     out/table_C.tex   LaTeX tabular
%     out/table_C.txt   plain-text mirror
%   Metric: robust fraction per (model, eps, tool, algorithm) + mean time.

    p = inputParser;
    addParameter(p, 'resultsDir', ...
        fullfile(toolcomparison_root(), 'mnist_resnet', 'results'));
    addParameter(p, 'outDir', fullfile(toolcomparison_root(), 'tables', 'out'));
    parse(p, varargin{:});
    opts = p.Results;
    if ~isfolder(opts.outDir), mkdir(opts.outDir); end

    u = tool_utils();
    models = {'mnist_resnet8'};

    header = {"Model","\\epsilon","Tool","Algorithm","V","X","?","E","Mean t (s)","Total t (s)","PAR-2 (s)"};
    rows = {};
    txtLines = strings(0,1);
    for i = 1:numel(models)
        model   = models{i};
        matFile = fullfile(opts.resultsDir, sprintf("expC_%s.mat", model));
        R = u.load_results(matFile);
        if isempty(R), continue; end
        epsVals = extract_eps(R.instance_id);
        R.eps = epsVals;

        keys = unique(R(:, {'eps','tool','algorithm'}));
        for k = 1:height(keys)
            e    = keys.eps(k);
            tool = keys.tool(k);
            alg  = keys.algorithm(k);
            sel  = R(R.eps==e & R.tool==tool & R.algorithm==alg & R.benchmark==string(model), :);
            n     = height(sel);
            v     = sum(sel.status == "verified");
            x     = sum(sel.status == "violated");
            unk   = sum(sel.status == "unknown");
            errN  = sum(sel.status == "error");
            tlist = sel.time(sel.status ~= "timeout" & sel.status ~= "error" & ~isnan(sel.time));
            mt    = mean(tlist);
            tt    = sum(tlist);
            % PAR-2 score: unsolved instances counted at 2 * timeout. Use the
            % canonical timeout from this group (uniform per row in this half).
            tEff  = median(sel.timeout(~isnan(sel.timeout)));
            if isnan(tEff), tEff = 0; end
            par2  = u.par2(sel.time, sel.status, tEff);
            rows{end+1} = { model, sprintf("%.4f", e), tool, alg, ...
                sprintf("%d", v), sprintf("%d", x), sprintf("%d", unk), sprintf("%d", errN), ...
                u.format_time(mt), u.format_time(tt), u.format_time(par2) }; %#ok<AGROW>
            txtLines(end+1,1) = sprintf("%-15s eps=%.4f %-14s %-22s V=%2d X=%2d ?=%2d E=%2d (n=%2d) mean=%s tot=%s PAR-2=%s", ...
                model, e, tool, alg, v, x, unk, errN, n, u.format_time(mt), u.format_time(tt), u.format_time(par2)); %#ok<AGROW>
        end
    end

    u.emit_latex_table( fullfile(opts.outDir,'table_C.tex'), ...
        header, rows, ...
        "ToolComparison ResNet head-to-head: NNV vs MathWorks AIVL (additionLayer, R2024b+).", ...
        "tab:toolcomparison-mnist-resnet");
    fid = fopen(fullfile(opts.outDir,'table_C.txt'),'w');
    fprintf(fid, "%s\n", strjoin(txtLines, newline));
    fclose(fid);
    fprintf("Wrote %s and %s\n", ...
        fullfile(opts.outDir,'table_C.tex'), ...
        fullfile(opts.outDir,'table_C.txt'));
end

function eps = extract_eps(ids)
    eps = nan(numel(ids),1);
    for i = 1:numel(ids)
        tok = regexp(ids(i), "eps=([0-9.]+)", "tokens", "once");
        if ~isempty(tok), eps(i) = str2double(tok{1}); end
    end
end

function r = toolcomparison_root()
    r = fileparts(fileparts(mfilename('fullpath')));
end
