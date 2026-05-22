function run_experiments(toolcomparisonMode)
%RUN_EXPERIMENTS  Drive the NNV3.0 suite in the current MATLAB session.
%
% Used by run_smoke.m and run_full.m for the nnv3.0-online flow, where
% bash run_all.sh cannot be used (matlab -batch can't see the browser
% sign-in licence).
%
%   run_experiments('smoke')   % ToolComparison ~12 min (NNV-only sanity)
%   run_experiments('full')    % ToolComparison ~3-5 h (Tables 5, 6, 7)
%
% After all experiments finish, results are consolidated to
% NNV3.0/repeatability_logs/ and a top-level PAPER_COMPARISON.md is
% written there listing every output file plus which paper table /
% figure it backs. This mirrors the consolidation bash run_all.sh
% performs on the network-licence path so a single host bind mount
% captures everything a reviewer needs.
%
% ProbVer is placed last because in this single-session flow an
% OOM-SIGKILL inside cp-star reachability would tear down the entire
% MATLAB session (and the cached licence) -- putting it last preserves
% the other experiments' results. ProbVer also auto-skips when no GPU
% is visible to the container (cp-star requires CUDA).

    NNV3_EXAMPLES = '/home/matlab/nnv/code/nnv/examples/NNV3.0';
    LOG_DIR       = fullfile(NNV3_EXAMPLES, 'repeatability_logs');
    FIG_DIR       = fullfile(LOG_DIR, 'figures');
    if ~isfolder(FIG_DIR), mkdir(FIG_DIR); end

    warning('off','backtrace');
    warning('off','nnet_cnn_onnx:onnx:WarnAPIDeprecation');
    warning('off','nnet_cnn_onnx:onnx:InputLabelMismatch');
    warning('off','nnet_cnn:internal:cnn:analyzer:NetworkAnalyzer:NetworkHasWarnings');
    warning('off','MATLAB:mpath:nameNonexistentOrNotADirectory');
    warning('off','MATLAB:linkaxes:RequireDataAxes');
    addpath(genpath('/home/matlab/nnv/code/nnv'));
    try, parallel.gpu.enableCUDAForwardCompatibility(true); catch; end

    % Suppress figure pop-ups while keeping saveas/exportgraphics/print
    % calls in the experiment scripts fully functional. Invisible figures
    % render and save normally; they just don't appear on the browser
    % MATLAB desktop. Restored on exit so the user's interactive session
    % is unchanged after the runner returns.
    prevFigVis = get(0, 'DefaultFigureVisible');
    restoreFigVis = onCleanup(@() set(0, 'DefaultFigureVisible', prevFigVis)); %#ok<NASGU>
    set(0, 'DefaultFigureVisible', 'off');

    setenv('TOOLCOMPARISON_MODE', toolcomparisonMode);

    hasGPU = false;
    try
        hasGPU = (gpuDeviceCount > 0);
    catch
    end

    experiments = { ...
        'fairnnv',        'FairNNV',        'run_fairnnv.m',           true;   ...
        'gnnv',           'GNNV',           'run_gnn_experiments.m',   true;   ...
        'videostar',      'VideoStar',      'run_zoomin_4f.m',         true;   ...
        'modelstar',      'ModelStar',      'run_expt_for_compute.m',  true;   ...
        'toolcomparison', 'ToolComparison', 'run_toolcomparison.m',    true;   ...
        'probver',        'ProbVer',        'run_probver.m',           hasGPU; ...
    };

    n = size(experiments,1);
    status = cell(n,1);
    wall   = zeros(n,1);

    for k = 1:n
        name = experiments{k,1}; subdir = experiments{k,2};
        entry = experiments{k,3}; enabled = experiments{k,4};
        if ~enabled
            fprintf(['\n=== %-15s SKIPPED (no GPU; cp-star reachability ' ...
                     'requires CUDA -- restart the container with --gpus all ' ...
                     'to enable) ===\n'], name);
            status{k} = 'skipped';
            wall(k) = 0;
            continue;
        end
        fprintf('\n=== %-15s start: %s ===\n', name, datestr(now,'yyyy-mm-ddTHH:MM:SS'));
        t0 = tic;
        try
            run_isolated(fullfile(NNV3_EXAMPLES, subdir), entry, FIG_DIR, name);
            status{k} = 'ok';
        catch ME
            status{k} = sprintf('failed: %s', ME.message);
            fprintf(2, '=== %-15s FAILED: %s\n', name, ME.message);
        end
        wall(k) = toc(t0);
        fprintf('=== %-15s %s in %.0fs ===\n', name, status{k}, wall(k));
    end

    consolidate_outputs(NNV3_EXAMPLES, LOG_DIR);
    write_paper_comparison_md(LOG_DIR, experiments, status, wall, toolcomparisonMode);

    fprintf('\n=========================== SUMMARY ============================\n');
    fprintf('%-16s %-10s %s\n','experiment','wall(s)','status');
    fprintf('%-16s %-10s %s\n','----------','-------','------');
    for k = 1:n
        fprintf('%-16s %-10.0f %s\n', experiments{k,1}, wall(k), status{k});
    end
    fprintf('================================================================\n');
    fprintf('Consolidated outputs:    %s/results/\n', LOG_DIR);
    fprintf('Paper comparison index:  %s/PAPER_COMPARISON.md\n', LOG_DIR);
    if strcmp(toolcomparisonMode, 'full')
        fprintf('Headline paper tables:   %s/results/ToolComparison/tables/table_main.{tex,txt}\n', LOG_DIR);
        fprintf('   (compare against `aivl_comparison.tex` / Tables 5, 6, 7 in the paper)\n');
    end
    fprintf('\nTo extract to the host (from a separate PowerShell):\n');
    fprintf('  docker cp nnv3-setup:%s .\\repeatability_logs\n', LOG_DIR);

    total_wall = sum(wall);
    fprintf('\n');
    fprintf('################################################################\n');
    fprintf('################################################################\n');
    fprintf('##                                                            ##\n');
    fprintf('##              NNV 3.0 %-4s RUN COMPLETED                    ##\n', upper(toolcomparisonMode));
    fprintf('##              Total wall: %5.0f s  (%5.1f min)              ##\n', total_wall, total_wall/60);
    fprintf('##                                                            ##\n');
    fprintf('################################################################\n');
    fprintf('################################################################\n');
    fprintf('\nThe run has finished. You can now:\n');
    fprintf('  1. Open the paper comparison index inside the container:\n');
    fprintf('       %s/PAPER_COMPARISON.md\n', LOG_DIR);
    fprintf('  2. Or extract the whole repeatability_logs tree to the host\n');
    fprintf('     with the docker cp command shown above.\n');
    fprintf('  3. Exit the browser MATLAB session / Ctrl-C the docker run\n');
    fprintf('     terminal when you are done.\n\n');

    print_paper_artefacts(LOG_DIR);

    cd(NNV3_EXAMPLES);
end

% -------------------------------------------------------------------------
function run_isolated(entryDir, entry, figDir, expLabel)
% Run a single experiment script in this function's local workspace so
% any variables it creates (notably `config`, which several runners
% guard with `if ~exist('config','var')` and then *inherit* any
% pre-existing value from) cannot leak into sibling experiments in the
% shared single-MATLAB-session online flow.
%
% Without this isolation, FairNNV sets config.modelsDir = .../FairNNV/models;
% the variable persists across experiments and VideoStar's
% run_zoomin_4f.m then resolves `zoomin_4f.onnx` against FairNNV/models/
% and fails with "Model not found". The bash run_all.sh path doesn't hit
% this because every experiment runs in a fresh matlab -batch process.
%
% After the script returns, any figures it left open are captured as
% PNG + PDF under figDir/ (DefaultFigureVisible='off' from the prelude
% keeps them off-screen, but exportgraphics still renders them) and
% then closed so they don't leak into the next experiment. Used to
% surface paper Fig. 2 (ModelStar) which is otherwise drawn but never
% saved by the per-experiment script.

    cd(entryDir);
    run(entry);

    if nargin >= 4 && ~isempty(figDir)
        figs = findall(groot, 'Type', 'figure');
        for fk = 1:numel(figs)
            base = fullfile(figDir, sprintf('%s_fig%d', expLabel, fk));
            try, exportgraphics(figs(fk), [base '.png'], 'Resolution', 200); catch; end
            try, exportgraphics(figs(fk), [base '.pdf']); catch; end
        end
        close all force
    end
end

% -------------------------------------------------------------------------
function print_paper_artefacts(logDir)
% Render the paper's Figs 2-3 and Tables 2-7 from the consolidated
% experiment outputs and tee the rendered output to both the MATLAB
% Command Window AND repeatability_logs/PAPER_TABLES.txt so reviewers
% have a single file to compare against the ATVA 2026 PDF.
%
% Paper -> artefact map:
%   Fig. 2  -> ModelStar single-layer weight-perturbation (paper Fig. 2)
%             Captured by run_isolated post-step from the figure ModelStar's
%             plot_results leaves open. The conv_expt_any_layer.m caller
%             passes reverse_order=0 so the columns appear as fc_6|fc_5|fc_4
%             matching the paper.
%   Fig. 3  -> GNNV IEEE-24 power flow (paper Fig. 3)
%             Generated here from gnn_results.csv (verified-% vs node-eps,
%             one line per architecture) -- GNNV's per-experiment script
%             doesn't save it.
%   Table 2 -> VolumeStar / VideoStar; aggregate the per-sample CSVs into
%             per-epsilon (Ver, Unk, Avg time) rows.
%   Table 3 -> ProbVer; extract property number from the vnnlib filename.
%   Table 4 -> FairNNV; combine counterfactual + individual + timing CSVs
%             into the paper's (Small / Medium) x (VF% / Time) layout.
%   Tables 5+6+7 -> ToolComparison; loaded straight from per-benchmark
%             .mat files via the canonical tool_utils schema, split into
%             FC VNNLIB (ACAS p3/p4/RL), CNN VNNLIB (OVAL21, Collins RUL),
%             and MNIST-ResNet-8 (per epsilon, decoded from instance_id).

    resultsDir   = fullfile(logDir, 'results');
    figDir       = fullfile(logDir, 'figures');
    paperTxt     = fullfile(logDir, 'PAPER_TABLES.txt');

    mode    = getenv('TOOLCOMPARISON_MODE');
    if isempty(mode), mode = 'unknown'; end
    isSmoke = strcmp(mode, 'smoke');

    % Generate Fig. 3 before listing the figures section.
    try
        generate_fig3(resultsDir, figDir);
    catch ME
        fprintf(2, '  [fig 3 generation] %s\n', ME.message);
    end

    tfid = fopen(paperTxt, 'w');
    if tfid > 0
        fids = [1, tfid];
        cleanup = onCleanup(@() fclose(tfid)); %#ok<NASGU>
    else
        fids = 1;
    end

    dprintf(fids, '\n');
    dprintf(fids, '================================================================\n');
    dprintf(fids, '  PAPER ARTEFACTS (compare against the NNV3 ATVA 2026 PDF)\n');
    dprintf(fids, '  Run mode: %s\n', mode);
    dprintf(fids, '================================================================\n');

    dprintf(fids, '\nFIGURES saved as PNG + PDF under:\n  %s/\n', figDir);
    dprintf(fids, '\n  Fig. 2 -- ModelStar single-layer weight-perturbation (paper Fig. 2)\n');
    list_fig_matches(fids, figDir, 'modelstar_*', ...
        '    Source: results/ModelStar/MNIST_MLP.mat (per-cell verdicts)');
    dprintf(fids, '    Layer order: fc_6 | fc_5 | fc_4  (paper-aligned via reverse_order=0)\n');
    dprintf(fids, '\n  Fig. 3 -- GNNV IEEE-24 power flow (paper Fig. 3)\n');
    list_fig_matches(fids, figDir, 'fig_3_gnnv*', ...
        '    Source: results/GNNV/gnn_results.csv');
    dprintf(fids, '    Plot: verified-%% vs node-epsilon, one line per arch (GCN / SAGE / GINE-Conv)\n');

    render_table_2(fids, resultsDir);
    render_table_3(fids, resultsDir);
    render_table_4(fids, resultsDir);
    render_tables_5_6_7(fids, resultsDir, isSmoke);

    dprintf(fids, '\n================================================================\n');
    if tfid > 0
        dprintf(fids, 'Rendered tables saved to: %s\n', paperTxt);
        dprintf(fids, 'docker cp the parent repeatability_logs/ dir to get every artefact.\n');
    end
    dprintf(fids, '================================================================\n\n');
end

% -------------------------------------------------------------------------
function dprintf(fids, varargin)
%DPRINTF  fprintf to a vector of file IDs (e.g. [1 fid] for stdout+file).
    for k = 1:numel(fids)
        fprintf(fids(k), varargin{:});
    end
end

% -------------------------------------------------------------------------
function list_fig_matches(fids, figDir, prefix, sourceNote)
    if ~isfolder(figDir)
        dprintf(fids, '    (no figure directory yet -- experiment may not have run)\n');
        if nargin >= 4, dprintf(fids, '%s\n', sourceNote); end
        return;
    end
    info = [dir(fullfile(figDir, [prefix '.png'])); ...
            dir(fullfile(figDir, [prefix '.pdf']))];
    if isempty(info)
        dprintf(fids, '    (not auto-generated by this experiment)\n');
        if nargin >= 4, dprintf(fids, '%s\n', sourceNote); end
    else
        for k = 1:numel(info)
            dprintf(fids, '    %s\n', fullfile(info(k).folder, info(k).name));
        end
    end
end

% -------------------------------------------------------------------------
function generate_fig3(resultsDir, figDir)
% Plot verified-% vs node-epsilon, one line per architecture, from
% GNNV's gnn_results.csv. GNNV doesn't save its figure; without this
% step there's no Fig. 3 artefact on disk.

    csv = fullfile(resultsDir, 'GNNV', 'gnn_results.csv');
    if ~isfile(csv), return; end
    T = readtable(csv);

    archs = unique(T.Architecture, 'stable');
    if isempty(archs), return; end

    fig = figure('Visible', 'off', 'Position', [100 100 720 480]);
    hold on;
    markers = {'o','s','d','^','v','>'};
    for k = 1:numel(archs)
        Tk = T(strcmp(T.Architecture, archs{k}), :);
        [eps_sorted, idx] = sort(Tk.Node_Epsilon);
        pct_sorted = Tk.Pct_Verified(idx);
        if numel(eps_sorted) >= 1
            semilogx(eps_sorted, pct_sorted, '-', ...
                'Marker', markers{mod(k-1, numel(markers))+1}, ...
                'LineWidth', 2, 'MarkerSize', 8, ...
                'DisplayName', upper(archs{k}));
        end
    end
    xlabel('Node \epsilon');
    ylabel('Verified voltage nodes (%)');
    title('Fig. 3 -- GNNV IEEE-24 power flow verification');
    legend('Location', 'best');
    grid on;
    ylim([0 105]);

    if ~isfolder(figDir), mkdir(figDir); end
    try, exportgraphics(fig, fullfile(figDir, 'fig_3_gnnv.png'), 'Resolution', 200); catch; end
    try, exportgraphics(fig, fullfile(figDir, 'fig_3_gnnv.pdf')); catch; end
    close(fig);
end

% -------------------------------------------------------------------------
function render_table_2(fids, resultsDir)
% Paper Table 2: VolumeStar verification on ZoomIn-4f.
% Aggregate per-sample eps=*.csv into per-epsilon (Ver, Unk, Avg time).
% Result codes: 1 = verified, 2 = unknown.

    dprintf(fids, '\n--- Table 2 -- VolumeStar / VideoStar (paper Tab. 2) ---\n');
    dprintf(fids, '10 samples per epsilon, 30-min timeout, relax algorithm\n\n');

    epsCases = {'1/255', '2/255', '3/255'};
    csvFiles = {'eps=1_255.csv', 'eps=2_255.csv', 'eps=3_255.csv'};

    dprintf(fids, '  %-9s | %-4s | %-4s | %-14s\n', 'epsilon', 'Ver.', 'Unk.', 'Avg. Time (s)');
    dprintf(fids, '  %s\n', '----------+------+------+---------------');
    for k = 1:numel(epsCases)
        f = fullfile(resultsDir, 'VideoStar', csvFiles{k});
        if ~isfile(f)
            dprintf(fids, '  %-9s | %-4s | %-4s | %-14s\n', epsCases{k}, '?', '?', '(missing)');
            continue;
        end
        try
            T = readtable(f);
            nV = sum(T.Result == 1);
            nU = sum(T.Result == 2);
            avgT = mean(T.Time);
            dprintf(fids, '  %-9s | %4d | %4d | %14.2f\n', epsCases{k}, nV, nU, avgT);
        catch ME
            dprintf(fids, '  %-9s | %-4s | %-4s | (parse error: %s)\n', epsCases{k}, '?', '?', ME.message);
        end
    end
end

% -------------------------------------------------------------------------
function render_table_3(fids, resultsDir)
% Paper Table 3: Probabilistic verification on TinyYOLO.
% Replace onnx/vnnlib columns with property number extracted from
% the vnnlib filename ("TinyYOLO_prop_000101_eps_1_255.vnnlib" -> Prop 101).

    dprintf(fids, '\n--- Table 3 -- Probabilistic verification (paper Tab. 3) ---\n');
    dprintf(fids, 'CP-Star reachability, surrogate 99.9%% coverage / 99.9%% confidence, GPU\n\n');

    f = fullfile(resultsDir, 'ProbVer', 'results_summary.csv');
    if ~isfile(f)
        dprintf(fids, '  (results_summary.csv not found at %s)\n', f);
        return;
    end
    try
        T = readtable(f);
    catch ME
        dprintf(fids, '  (could not read: %s)\n', ME.message);
        return;
    end

    dprintf(fids, '  %-10s | %-7s | %-9s | %-7s\n', 'Property', 'epsilon', 'Time (s)', 'Result');
    dprintf(fids, '  %s\n', '-----------+---------+-----------+--------');
    for k = 1:height(T)
        v = char(T.vnnlib{k});
        propTok = regexp(v, 'prop_0*(\d+)', 'tokens', 'once');
        epsTok  = regexp(v, 'eps_(\d+)_(\d+)', 'tokens', 'once');
        if isempty(propTok), propStr = ['#' num2str(T.index(k))]; else, propStr = ['Prop ' propTok{1}]; end
        if isempty(epsTok),  epsStr  = '?'; else, epsStr = [epsTok{1} '/' epsTok{2}]; end
        statusStr = upper(char(T.status{k}));
        dprintf(fids, '  %-10s | %-7s | %9.2f | %-7s\n', propStr, epsStr, T.time(k), statusStr);
    end
end

% -------------------------------------------------------------------------
function render_table_4(fids, resultsDir)
% Paper Table 4: FairNNV on Adult Census.
% Combine counterfactual + individual fairness CSVs with the timing CSV
% into the paper's two-row-per-model (VF % / Time s) layout. Use only
% the latest timestamped run (the consolidator may carry duplicates).
% Display labels: AC-1 -> Small, AC-3 -> Medium.

    dprintf(fids, '\n--- Table 4 -- FairNNV on Adult Census (paper Tab. 4) ---\n');
    dprintf(fids, 'CF = counterfactual fairness (epsilon=0); IF = individual fairness\n\n');

    fnDir = fullfile(resultsDir, 'FairNNV');
    cfPath  = latest_match(fnDir, 'counterfactual_*.csv');
    indPath = latest_match(fnDir, 'individual_*.csv');
    timPath = latest_match(fnDir, 'timing_*.csv');

    if isempty(cfPath) || isempty(indPath)
        dprintf(fids, '  (FairNNV CSVs not found under %s)\n', fnDir);
        return;
    end

    try
        cfT  = readtable(cfPath);
        indT = readtable(indPath);
    catch ME
        dprintf(fids, '  (could not read CSVs: %s)\n', ME.message);
        return;
    end
    timT = [];
    if ~isempty(timPath), try, timT = readtable(timPath); catch; end; end

    eps_list = unique(indT.Epsilon, 'stable');
    if isempty(eps_list)
        dprintf(fids, '  (no epsilons in individual_*.csv)\n');
        return;
    end

    % Header
    head = sprintf('  %-9s %-10s   %-8s', 'Model', 'Metric', 'CF(e=0)');
    for k = 1:numel(eps_list)
        head = [head sprintf(' %-8s', sprintf('IF(%.2f)', eps_list(k)))]; %#ok<AGROW>
    end
    dprintf(fids, '%s\n', head);
    dprintf(fids, '  %s\n', repmat('-', 1, length(head)-2));

    modelMap = {{'AC-1','Small'}, {'AC-3','Medium'}};
    for mi = 1:numel(modelMap)
        key  = modelMap{mi}{1};
        disp_= modelMap{mi}{2};

        cfRow = cfT(strcmp(cfT.Model, key), :);
        cfVF  = NaN; if ~isempty(cfRow), cfVF = cfRow.FairPercent(1); end

        % VF row
        line = sprintf('  %-9s %-10s   %8.0f', disp_, 'VF (%)', cfVF);
        for k = 1:numel(eps_list)
            r = indT(strcmp(indT.Model, key) & indT.Epsilon == eps_list(k), :);
            v = NaN; if ~isempty(r), v = r.FairPercent(1); end
            line = [line sprintf(' %8.0f', v)]; %#ok<AGROW>
        end
        dprintf(fids, '%s\n', line);

        % Time row (only if timing CSV present)
        if ~isempty(timT)
            cfTimeRow = timT(strcmp(timT.Model, key) & timT.Epsilon == 0, :);
            cfTime = NaN;
            if ~isempty(cfTimeRow) && ismember('AvgTimePerSample', timT.Properties.VariableNames)
                cfTime = cfTimeRow.AvgTimePerSample(1);
            end
            tline = sprintf('  %-9s %-10s   %8.2f', '', 'Time (s)', cfTime);
            for k = 1:numel(eps_list)
                r = timT(strcmp(timT.Model, key) & timT.Epsilon == eps_list(k), :);
                v = NaN;
                if ~isempty(r) && ismember('AvgTimePerSample', timT.Properties.VariableNames)
                    v = r.AvgTimePerSample(1);
                end
                tline = [tline sprintf(' %8.2f', v)]; %#ok<AGROW>
            end
            dprintf(fids, '%s\n', tline);
        end
    end

    dprintf(fids, '\n  Note: ''Small'' aggregates AC-1 rows, ''Medium'' aggregates AC-3 rows.\n');
    if isempty(timT)
        dprintf(fids, '  Note: per-sample timing data not found (timing_*.csv missing).\n');
    end
end

% -------------------------------------------------------------------------
function render_tables_5_6_7(fids, resultsDir, isSmoke)
% Paper Tables 5, 6, 7: ToolComparison NNV vs AIVL.
% Load per-benchmark .mat files via the canonical tool_utils schema
% (tool/benchmark/instance_id/status/time/algorithm/timeout/note),
% aggregate, and render the paper's three separate tables. For MNIST-
% ResNet-8 (Table 7), decode epsilon from the instance_id string
% ("img%d|eps=%s") to get per-epsilon mean times.

    smokeTag = '';
    if isSmoke, smokeTag = ' (SMOKE -- ~1 instance per cell)'; end

    % --- Table 5: FC VNNLIB ---
    dprintf(fids, '\n--- Table 5 -- Tool comparison on FC VNNLIB benchmarks (paper Tab. 5)%s ---\n', smokeTag);
    dprintf(fids, 'Per cell: V verified, X violated, U unknown, T timeout, S mean time on solved (s)\n');
    render_tc_grouped(fids, resultsDir, ...
        {'acas_xu_p3', 'acas_xu_p4', 'rl'}, {'ACAS p3', 'ACAS p4', 'RL'});

    % --- Table 6: CNN VNNLIB ---
    dprintf(fids, '\n--- Table 6 -- Tool comparison on CNN VNNLIB benchmarks (paper Tab. 6)%s ---\n', smokeTag);
    dprintf(fids, 'Per cell: V verified, X violated, U unknown, T timeout, S mean time on solved (s)\n');
    render_tc_grouped(fids, resultsDir, ...
        {'oval21', 'collins_rul'}, {'OVAL21', 'Collins RUL'});

    % --- Table 7: MNIST-ResNet-8 per epsilon ---
    dprintf(fids, '\n--- Table 7 -- MNIST-ResNet-8 robustness per epsilon (paper Tab. 7)%s ---\n', smokeTag);
    render_table_7_resnet(fids, resultsDir);

    if isSmoke
        dprintf(fids, '\n  Smoke mode is a structural sanity check, not a paper-faithful\n');
        dprintf(fids, '  reproduction. For the paper numbers, run run_full.m\n');
        dprintf(fids, '  (~3-5 h ToolComparison wall).\n');
    end
end

% -------------------------------------------------------------------------
function render_tc_grouped(fids, resultsDir, benches, labels)
% Aggregate each benchmark's .mat by (tool, algorithm) and print a
% side-by-side table with V/X/U/T/S columns per benchmark.

    aggregated = cell(1, numel(benches));
    Ns = zeros(1, numel(benches));
    for b = 1:numel(benches)
        matFile = fullfile(resultsDir, 'ToolComparison', [benches{b} '.mat']);
        R = load_tc_results(matFile);
        if ~istable(R) || isempty(R)
            aggregated{b} = []; Ns(b) = 0; continue;
        end
        aggregated{b} = aggregate_by_algo(R);
        Ns(b) = max(arrayfun(@(k) sum(R.tool == aggregated{b}.tool(k) & ...
                                       R.algorithm == aggregated{b}.algorithm(k)), ...
                              1:height(aggregated{b})));
    end

    if all(cellfun(@isempty, aggregated))
        dprintf(fids, '  (no ToolComparison results found under %s)\n', ...
            fullfile(resultsDir, 'ToolComparison'));
        return;
    end

    % Header rows
    cell_w = 26;
    headRow1 = sprintf('  %-26s', 'Algorithm');
    headRow2 = sprintf('  %-26s', '');
    for b = 1:numel(benches)
        headRow1 = [headRow1 ' | ' pad_center(sprintf('%s (N=%d)', labels{b}, Ns(b)), cell_w)]; %#ok<AGROW>
        headRow2 = [headRow2 ' | ' pad_center(' V   X   U   T      S', cell_w)]; %#ok<AGROW>
    end
    dprintf(fids, '\n%s\n%s\n', headRow1, headRow2);
    dprintf(fids, '  %s\n', repmat('-', 1, length(headRow1)-2));

    % Union of (tool, algorithm) keys across all benchmarks, AIVL first.
    keys = strings(0,1);
    for b = 1:numel(benches)
        if isempty(aggregated{b}), continue; end
        for r = 1:height(aggregated{b})
            keys(end+1, 1) = aggregated{b}.tool(r) + "/" + aggregated{b}.algorithm(r); %#ok<AGROW>
        end
    end
    keys = unique(keys, 'stable');
    aivl_k = keys(startsWith(keys, "aivl/"));
    nnv_k  = sort(keys(startsWith(keys, "nnv/")));
    keys   = [aivl_k; nnv_k];

    for ki = 1:numel(keys)
        parts = strsplit(char(keys(ki)), '/');
        toolS = parts{1}; algoS = parts{2};
        rowLabel = sprintf('%-4s %s', toolS, algoS);
        line = sprintf('  %-26s', rowLabel);
        for b = 1:numel(benches)
            agg = aggregated{b};
            if isempty(agg)
                line = [line ' | ' pad_center('--', cell_w)]; %#ok<AGROW>
                continue;
            end
            r = agg(agg.tool == toolS & agg.algorithm == algoS, :);
            if isempty(r)
                cell_s = '--';
            else
                if isnan(r.MeanTime(1)), ts = '  --'; else, ts = sprintf('%6.2f', r.MeanTime(1)); end
                cell_s = sprintf('%2d  %2d  %2d  %2d  %s', r.V(1), r.X(1), r.U(1), r.T(1), ts);
            end
            line = [line ' | ' pad_center(cell_s, cell_w)]; %#ok<AGROW>
        end
        dprintf(fids, '%s\n', line);
    end
end

% -------------------------------------------------------------------------
function render_table_7_resnet(fids, resultsDir)
    matFile = fullfile(resultsDir, 'ToolComparison', 'mnist_resnet8.mat');
    R = load_tc_results(matFile);
    if ~istable(R) || isempty(R)
        dprintf(fids, '  (mnist_resnet8.mat not found under %s)\n', ...
            fullfile(resultsDir, 'ToolComparison'));
        return;
    end

    % Decode epsilon from instance_id = "img%d|eps=%s"
    eps_vec = nan(height(R), 1);
    for k = 1:height(R)
        tok = regexp(char(R.instance_id(k)), 'eps=([\d.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), eps_vec(k) = str2double(tok{1}); end
    end
    R.eps = eps_vec;
    R = R(~isnan(R.eps), :);
    if isempty(R)
        dprintf(fids, '  (no eps-tagged instance_ids in mnist_resnet8.mat)\n');
        return;
    end

    eps_unique = sort(unique(R.eps));
    eps_labels = arrayfun(@(e) format_eps(e), eps_unique, 'UniformOutput', false);

    % Determine tool/algorithm pairs in canonical order: AIVL first.
    pairs = unique(R(:, {'tool','algorithm'}), 'rows', 'stable');
    aivlMask = strcmp(string(pairs.tool), 'aivl');
    pairs = [pairs(aivlMask, :); pairs(~aivlMask, :)];

    dprintf(fids, '%d images per epsilon (full mode renders all eps; smoke mode is single-eps).\n\n', ...
        round(height(R) / max(1, numel(eps_unique) * height(pairs))));

    head = sprintf('  %-30s', 'Tool / Algorithm');
    for k = 1:numel(eps_unique), head = [head sprintf(' | %8s', eps_labels{k})]; end %#ok<AGROW>
    head2 = sprintf('  %-30s', '');
    for k = 1:numel(eps_unique), head2 = [head2 sprintf(' | %8s', 'mean t(s)')]; end %#ok<AGROW>
    dprintf(fids, '%s\n%s\n', head, head2);
    dprintf(fids, '  %s\n', repmat('-', 1, length(head)-2));

    for pi = 1:height(pairs)
        toolS = string(pairs.tool(pi)); algoS = string(pairs.algorithm(pi));
        line = sprintf('  %-30s', sprintf('%-4s %s', toolS, algoS));
        for ei = 1:numel(eps_unique)
            e = eps_unique(ei);
            sub = R(R.tool == toolS & R.algorithm == algoS & abs(R.eps - e) < 1e-9, :);
            solved = sub.status == "verified" | sub.status == "violated";
            if any(solved), mt = mean(sub.time(solved)); else, mt = NaN; end
            if isnan(mt), line = [line sprintf(' | %8s', '--')]; %#ok<AGROW>
            else,         line = [line sprintf(' | %8.2f', mt)]; %#ok<AGROW>
            end
        end
        dprintf(fids, '%s\n', line);
    end
end

% -------------------------------------------------------------------------
function agg = aggregate_by_algo(R)
    keys = unique(R(:, {'tool','algorithm'}), 'rows', 'stable');
    n = height(keys);
    agg = table('Size', [n 7], ...
        'VariableTypes', {'string','string','double','double','double','double','double'}, ...
        'VariableNames', {'tool','algorithm','V','X','U','T','MeanTime'});
    for k = 1:n
        toolS = string(keys.tool(k));
        algoS = string(keys.algorithm(k));
        sub   = R(string(R.tool) == toolS & string(R.algorithm) == algoS, :);
        agg.tool(k) = toolS; agg.algorithm(k) = algoS;
        agg.V(k) = sum(string(sub.status) == "verified");
        agg.X(k) = sum(string(sub.status) == "violated");
        agg.U(k) = sum(string(sub.status) == "unknown");
        agg.T(k) = sum(string(sub.status) == "timeout");
        solved   = string(sub.status) == "verified" | string(sub.status) == "violated";
        if any(solved), agg.MeanTime(k) = mean(sub.time(solved));
        else,           agg.MeanTime(k) = NaN;
        end
    end
end

% -------------------------------------------------------------------------
function R = load_tc_results(matFile)
    R = [];
    if ~isfile(matFile), return; end
    try
        S = load(matFile);
        f = fieldnames(S);
        if isempty(f), return; end
        R = S.(f{1});
    catch
        R = [];
    end
end

% -------------------------------------------------------------------------
function p = latest_match(d, pat)
    p = '';
    if ~isfolder(d), return; end
    info = dir(fullfile(d, pat));
    if isempty(info), return; end
    [~, idx] = max([info.datenum]);
    p = fullfile(info(idx).folder, info(idx).name);
end

% -------------------------------------------------------------------------
function s = format_eps(e)
    if abs(e * 255 - round(e*255)) < 1e-3
        s = sprintf('%d/255', round(e*255));
    else
        s = sprintf('%.5g', e);
    end
end

% -------------------------------------------------------------------------
function s = pad_center(str, w)
    str = char(str);
    if length(str) >= w, s = str(1:w); return; end
    pad = w - length(str);
    left = floor(pad/2); right = pad - left;
    s = [repmat(' ', 1, left) str repmat(' ', 1, right)];
end

% -------------------------------------------------------------------------
function consolidate_outputs(examplesDir, logDir)
% Mirror bash run_all.sh's consolidation: copy each experiment's latest
% timestamped results subdir + the ToolComparison results and rendered
% tables into a single tree at logDir/results/.

    if ~isfolder(logDir), mkdir(logDir); end
    resultsDir = fullfile(logDir, 'results');
    if ~isfolder(resultsDir), mkdir(resultsDir); end

    expDirs = {'FairNNV', 'ProbVer', 'GNNV', 'VideoStar', 'ModelStar'};
    for k = 1:numel(expDirs)
        src = fullfile(examplesDir, expDirs{k}, 'results');
        if ~isfolder(src), continue; end
        d = dir(src);
        sub = d([d.isdir] & ~ismember({d.name}, {'.','..'}));
        if ~isempty(sub)
            [~, idx] = max([sub.datenum]);
            srcDir = fullfile(src, sub(idx).name);
        else
            srcDir = src;
        end
        dst = fullfile(resultsDir, expDirs{k});
        try, copyfile(srcDir, dst, 'f'); catch ME, ...
            fprintf(2,'[consolidate] %s skipped: %s\n', expDirs{k}, ME.message); end
    end

    tcResults = fullfile(examplesDir, 'ToolComparison', 'results');
    if isfolder(tcResults)
        tcDst = fullfile(resultsDir, 'ToolComparison');
        if ~isfolder(tcDst), mkdir(tcDst); end
        try, copyfile(fullfile(tcResults, '*'), tcDst, 'f'); catch; end
    end
    tcTables = fullfile(examplesDir, 'ToolComparison', 'tables', 'out');
    if isfolder(tcTables)
        tcTablesDst = fullfile(resultsDir, 'ToolComparison', 'tables');
        if ~isfolder(tcTablesDst), mkdir(tcTablesDst); end
        try, copyfile(fullfile(tcTables, '*'), tcTablesDst, 'f'); catch; end
    end
end

% -------------------------------------------------------------------------
function write_paper_comparison_md(logDir, experiments, status, wall, mode)
% Write a single index file (PAPER_COMPARISON.md) at the consolidated
% results root, mapping each experiment to its primary output and the
% paper artefact it backs. Reviewers open one file to find everything.

    if ~isfolder(logDir), mkdir(logDir); end
    mdPath = fullfile(logDir, 'PAPER_COMPARISON.md');
    fid = fopen(mdPath, 'w');
    if fid < 0, return; end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    primary = {
        'fairnnv',        'results/FairNNV/timing_*.csv + counterfactual/individual reports', 'paper §FairNNV (fairness fractions table)';
        'gnnv',           'results/GNNV/gnn_results.csv + dashboard.{png,pdf}',               'paper §GNNV (per-architecture voltage bounds)';
        'videostar',      'results/VideoStar/eps=*.csv',                                       'paper §VideoStar (ZoomIn-4f epsilon sweep)';
        'modelstar',      'results/ModelStar/MNIST_MLP.mat + heatmap',                         'paper §ModelStar (weight-perturbation heatmap)';
        'toolcomparison', 'results/ToolComparison/tables/table_main.{tex,txt}',                'paper Tables 5, 6, 7 (NNV-vs-AIVL comparison)';
        'probver',        'results/ProbVer/results_summary.csv',                               'paper §ProbVer (TinyYOLO cp-star verdicts)';
    };

    fprintf(fid, '# NNV 3.0 ATVA 2026 -- Paper Comparison Reference\n\n');
    fprintf(fid, 'Generated %s by `run_%s.m` (online-licence flow).\n', ...
        datestr(now,'yyyy-mm-ddTHH:MM:SS'), mode);
    fprintf(fid, 'Mode: **%s**\n\n', mode);
    fprintf(fid, 'All paths in this file are relative to the directory that contains\n');
    fprintf(fid, 'this `PAPER_COMPARISON.md` (i.e. `repeatability_logs/` inside the\n');
    fprintf(fid, 'container, or whatever you copied it to on the host).\n\n');

    fprintf(fid, '## Headline result (paper Tables 5, 6, 7)\n\n');
    fprintf(fid, '- LaTeX: `results/ToolComparison/tables/table_main.tex`\n');
    fprintf(fid, '- Text:  `results/ToolComparison/tables/table_main.txt`\n\n');
    if strcmp(mode, 'full')
        fprintf(fid, 'These are the NNV-vs-AIVL comparison the paper renders as\n');
        fprintf(fid, '`aivl_comparison.tex`. Compare row-by-row against the published\n');
        fprintf(fid, 'Tables 5 (ACAS Xu), 6 (RL / OVAL21 / Collins RUL), and 7\n');
        fprintf(fid, '(MNIST-ResNet-8). PAR-2 and per-row Timeout columns are reported.\n\n');
    else
        fprintf(fid, 'In **smoke** mode ToolComparison runs 1 instance per benchmark\n');
        fprintf(fid, 'for sanity-checking only; the rendered table is a structurally\n');
        fprintf(fid, 'correct but quantitatively reduced version of the paper''s\n');
        fprintf(fid, 'headline tables. Run `run_full.m` for the ~3-5 h reproduction.\n\n');
    end

    fprintf(fid, '## Per-experiment outputs\n\n');
    fprintf(fid, '| Experiment | Status | Wall (s) | Primary output(s) | Backs |\n');
    fprintf(fid, '|---|---|---:|---|---|\n');
    for k = 1:size(experiments,1)
        name = experiments{k,1};
        idx = find(strcmp(primary(:,1), name), 1);
        if isempty(idx)
            prim = '(see experiment subdir)'; backs = '-';
        else
            prim = primary{idx,2}; backs = primary{idx,3};
        end
        fprintf(fid, '| %s | %s | %.0f | `%s` | %s |\n', ...
            name, status{k}, wall(k), prim, backs);
    end
    fprintf(fid, '\n');

    fprintf(fid, '## Reproduction notes\n\n');
    fprintf(fid, '- Runs were performed in a single MATLAB browser session because\n');
    fprintf(fid, '  the online (signed-in) licence used by `nnv3.0-online` cannot be\n');
    fprintf(fid, '  consumed by `matlab -batch` (which `bash run_all.sh` invokes).\n');
    fprintf(fid, '- See `code/nnv/examples/NNV3.0/ONLINE_LICENSE.md` for the full\n');
    fprintf(fid, '  setup-and-run protocol used to produce these outputs.\n');
    fprintf(fid, '- For the network-licence (port@host) reproduction path, see the\n');
    fprintf(fid, '  root `README.md` and `code/nnv/examples/NNV3.0/README.md`.\n');
end
