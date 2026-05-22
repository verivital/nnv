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
% Surface the paper figures and tables at the end of the run so a
% reviewer can directly compare against the ATVA 2026 PDF without
% hunting through subdirectories.
%
% Paper -> artefact map (per the NNV3 ATVA 2026 paper):
%   Fig. 2  -> ModelStar single-layer weight-perturbation
%   Fig. 3  -> GNNV IEEE-24 power flow (3 architectures)
%   Table 2 -> VolumeStar / VideoStar verdicts (eps=*.csv)
%   Table 3 -> ProbVer TinyYOLO verdicts (results_summary.csv)
%   Table 4 -> FairNNV Adult Census (counterfactual + individual CSVs)
%   Table 5 -> ToolComparison FC VNNLIB  (table_A.txt; full mode)
%   Table 6 -> ToolComparison CNN VNNLIB (table_A.txt; full mode)
%   Table 7 -> ToolComparison MNIST-ResNet-8 (table_C.txt; full mode)
%   In smoke mode the full per-table renders aren't produced; we fall
%   back to ToolComparison's table_main.{tex,txt}.

    resultsDir = fullfile(logDir, 'results');
    figDir     = fullfile(logDir, 'figures');

    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('  PAPER ARTEFACTS (compare against the NNV3 ATVA 2026 PDF)\n');
    fprintf('================================================================\n');

    fprintf('\nFIGURES saved as PNG + PDF under:\n  %s/\n', figDir);
    fprintf('  Fig. 2 -- ModelStar single-layer weight-perturbation\n');
    list_fig_matches(figDir, 'modelstar_*', '    Source: results/ModelStar/MNIST_MLP.mat (per-cell verdicts)');
    fprintf('  Fig. 3 -- GNNV IEEE-24 power flow (GCN / SAGE / GINE-Conv)\n');
    list_fig_matches(figDir, 'gnnv_*', '    Source: results/GNNV/gnn_results.csv (verified percentages per arch x epsilon)');

    print_table_file('Table 2 -- VolumeStar / VideoStar (paper Sec. on VolumeStar)', ...
        fullfile(resultsDir, 'VideoStar', 'eps=1_255.csv'), ...
        fullfile(resultsDir, 'VideoStar', 'eps=2_255.csv'), ...
        fullfile(resultsDir, 'VideoStar', 'eps=3_255.csv'));

    print_table_file('Table 3 -- Probabilistic verification (paper Sec. on ProbVer)', ...
        fullfile(resultsDir, 'ProbVer', 'results_summary.csv'));

    print_table_glob('Table 4 -- FairNNV on Adult Census', ...
        fullfile(resultsDir, 'FairNNV'), {'counterfactual_*.csv', 'individual_*.csv'});

    tc_tables  = fullfile(resultsDir, 'ToolComparison', 'tables');
    table_A    = fullfile(tc_tables, 'table_A.txt');
    table_C    = fullfile(tc_tables, 'table_C.txt');
    table_main = fullfile(tc_tables, 'table_main.txt');
    if exist(table_A, 'file') || exist(table_C, 'file')
        print_table_file('Tables 5 + 6 -- Tool comparison (FC + CNN VNNLIB)', table_A);
        print_table_file('Table 7 -- MNIST-ResNet-8 robustness', table_C);
    elseif exist(table_main, 'file')
        print_table_file('Tables 5, 6, 7 -- Tool comparison (smoke mode -- table_main aggregates all three)', table_main);
    else
        fprintf('\n--- Tables 5, 6, 7 ---\n  (no ToolComparison table found at %s)\n', tc_tables);
    end

    fprintf('\n================================================================\n\n');
end

% -------------------------------------------------------------------------
function list_fig_matches(figDir, prefix, sourceNote)
    if ~isfolder(figDir)
        fprintf('    (no figure directory yet -- experiment may not have run)\n');
        if nargin >= 3, fprintf('%s\n', sourceNote); end
        return;
    end
    info = [dir(fullfile(figDir, [prefix '.png'])); dir(fullfile(figDir, [prefix '.pdf']))];
    if isempty(info)
        fprintf('    (not auto-generated by this experiment)\n');
        if nargin >= 3, fprintf('%s\n', sourceNote); end
    else
        for k = 1:numel(info)
            fprintf('    %s\n', fullfile(info(k).folder, info(k).name));
        end
    end
end

% -------------------------------------------------------------------------
function print_table_file(label, varargin)
    fprintf('\n--- %s ---\n', label);
    for k = 1:numel(varargin)
        f = varargin{k};
        if isempty(f), continue; end
        fprintf('\n  [%s]\n', f);
        if exist(f, 'file')
            try
                txt = fileread(f);
                fprintf('%s\n', txt);
            catch ME
                fprintf('    (could not read: %s)\n', ME.message);
            end
        else
            fprintf('    (not found)\n');
        end
    end
end

% -------------------------------------------------------------------------
function print_table_glob(label, baseDir, patterns)
    fprintf('\n--- %s ---\n', label);
    if ~isfolder(baseDir)
        fprintf('  (directory not found: %s)\n', baseDir);
        return;
    end
    for k = 1:numel(patterns)
        info = dir(fullfile(baseDir, patterns{k}));
        for j = 1:numel(info)
            f = fullfile(info(j).folder, info(j).name);
            fprintf('\n  [%s]\n', f);
            try
                txt = fileread(f);
                fprintf('%s\n', txt);
            catch
                fprintf('    (could not read)\n');
            end
        end
    end
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
