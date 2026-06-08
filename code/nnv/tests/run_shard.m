function run_shard(shardIdx, numShards, mode)
% RUN_SHARD  Run a deterministic 1-of-N shard of the NNV test suite (CI fan-out).
%
%   run_shard(k, N)             % standard shard k of N (see exclusions below)
%   run_shard(k, N, 'standard') % same as above
%   run_shard(k, N, 'highmem')  % ONLY the anticipated-heavy tests (big-RAM/self-hosted/ACCRE)
%   run_shard(k, N, 'all')      % everything CI can run (excludes only GPU + broken-under-runtests)
%
% Designed for GitHub Actions `strategy.matrix` via matlab-actions/run-command:
%   cd('code/nnv'); install; addpath(genpath(fullfile(pwd,'tests')));
%   run_shard(${{ matrix.shard }}, ${{ N }}, 'standard');
%
% Writes a JUnit XML report to code/nnv/test-results-ci/ (uploaded as a CI artifact)
% and calls assertSuccess(results) so the job FAILS if any test in the shard fails.
%
% PHILOSOPHY (maximize standard coverage): exclude from the standard ubuntu-latest
% shards ONLY what genuinely cannot run there. The "anticipated-heavy" tests are
% INCLUDED by default and only moved to confirmedOomPatterns once a real CI run proves
% they OOM the ~16 GB runner.

    if nargin < 3 || isempty(mode), mode = 'standard'; end
    if ischar(shardIdx) || isstring(shardIdx), shardIdx = str2double(shardIdx); end
    if ischar(numShards) || isstring(numShards), numShards = str2double(numShards); end
    assert(numShards >= 1 && shardIdx >= 1 && shardIdx <= numShards, ...
        'run_shard:args', 'Require 1 <= shardIdx (%g) <= numShards (%g).', shardIdx, numShards);

    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.XMLPlugin
    import matlab.unittest.Verbosity

    testsDir = fileparts(mfilename('fullpath'));

    % (a) GPU tests — no GPU on GitHub-hosted runners (error / incomplete there).
    gpuPatterns = {'_gpu_', 'toGPU'};

    % (b) Known broken UNDER runtests — these tests share variables across %% sections,
    %     but runtests runs each %% in a fresh workspace (they pass when run as a script,
    %     fail as unit tests: "Unrecognized function or variable 'results'/'N_SAMPLES'").
    %     Excluded until the sections are made self-contained. TODO(SLM soundness work).
    brokenUnderRuntests = {'test_SLM_layers_soundness', 'test_SoftmaxLayer_reach_soundness'};

    % (c) CONFIRMED to OOM the ~16 GB ubuntu-latest runner. Populate ONLY from a real CI
    %     run that proves OOM; keep minimal to maximize standard coverage. (Empty so far.)
    confirmedOomPatterns = { };

    % Anticipated-heavy (memory/time) tests — NOT excluded from standard; we run them in
    % CI and move only proven-OOM ones into confirmedOomPatterns. mode='highmem' runs just
    % these (for a targeted big-RAM/self-hosted runner or ACCRE SLURM array).
    anticipatedHeavyPatterns = {'VolumeStar', 'Conv3D', 'AveragePooling3D', ...
        'Segmentation', 'SEGNET', 'PixelClassification', 'reach_exact', 'comprehensive_network'};

    ciExclude = [gpuPatterns, brokenUnderRuntests, confirmedOomPatterns];

    suite = TestSuite.fromFolder(testsDir, 'IncludingSubfolders', true);
    names = string({suite.Name});

    switch lower(mode)
        case 'standard', sel = suite(~matchAny(names, ciExclude));
        case 'highmem',  sel = suite(matchAny(names, anticipatedHeavyPatterns) & ~matchAny(names, gpuPatterns));
        case 'all',      sel = suite(~matchAny(names, [gpuPatterns, brokenUnderRuntests]));
        otherwise, error('run_shard:mode', "mode must be 'standard', 'highmem', or 'all'.");
    end

    % Deterministic round-robin partition (interleaves heavy/light tests for balance).
    sel = sel(shardIdx:numShards:end);

    fprintf('== run_shard %d/%d mode=%s: running %d of %d total (excluded: %d gpu, %d broken, %d confirmed-OOM) ==\n', ...
        shardIdx, numShards, mode, numel(sel), numel(suite), ...
        nnz(matchAny(names, gpuPatterns)), nnz(matchAny(names, brokenUnderRuntests)), nnz(matchAny(names, confirmedOomPatterns)));

    outDir = fullfile(fileparts(testsDir), 'test-results-ci');   % code/nnv/test-results-ci
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    xmlFile = fullfile(outDir, sprintf('results-%s-shard%02d-of%02d.xml', lower(mode), shardIdx, numShards));
    fprintf('JUnit report -> %s\n', xmlFile);

    if isempty(sel)
        fprintf('No tests selected for this shard; writing empty report and exiting success.\n');
        fid = fopen(xmlFile, 'w'); fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n<testsuites/>\n'); fclose(fid);
        return;
    end

    runner = TestRunner.withTextOutput('OutputDetail', Verbosity.Concise);
    runner.addPlugin(XMLPlugin.producingJUnitFormat(xmlFile));
    results = runner.run(sel);

    fprintf('== shard %d/%d done: %d passed, %d failed, %d incomplete ==\n', ...
        shardIdx, numShards, nnz([results.Passed]), nnz([results.Failed]), nnz([results.Incomplete]));

    assertSuccess(results);
end

function tf = matchAny(names, patterns)
    tf = false(1, numel(names));
    for p = 1:numel(patterns)
        tf = tf | contains(names, patterns{p}, 'IgnoreCase', true);
    end
end
