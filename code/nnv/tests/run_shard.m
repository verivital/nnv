function run_shard(shardIdx, numShards, mode)
% RUN_SHARD  Run a deterministic 1-of-N shard of the NNV test suite (CI fan-out).
%
%   run_shard(k, N)             % standard shard k of N (EXCLUDES high-memory/GPU tests)
%   run_shard(k, N, 'standard') % same as above
%   run_shard(1, 1, 'highmem')  % ONLY the high-memory tests (run on a big-RAM/self-hosted runner or ACCRE)
%   run_shard(k, N, 'all')      % everything, no high-mem split (sharded)
%
% Designed for GitHub Actions `strategy.matrix` via matlab-actions/run-command:
%   cd('code/nnv'); install; addpath(genpath(fullfile(pwd,'tests')));
%   run_shard(${{ matrix.shard }}, ${{ N }}, 'standard');
%
% Writes a JUnit XML report to code/nnv/test-results-ci/ (uploaded as a CI artifact)
% and calls assertSuccess(results) so the job FAILS if any test in the shard fails.
%
% The standard ubuntu-latest GitHub-hosted runner has ~16 GB RAM (public repos), which
% several NNV tests exceed. Those are matched by highMemPatterns below, EXCLUDED from
% the standard shards, and run separately in 'highmem' mode on a capable runner.
% EDIT highMemPatterns as the full-suite duration/OOM report identifies new heavy tests.

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

    % --- High-memory / GPU test name patterns (case-insensitive substring on test Name) ---
    highMemPatterns = { ...
        'VolumeStar', ...                              % 3D/volumetric reachability (set/volume_star, soundness)
        'Conv3D', 'AveragePooling3D', ...              % 3D conv/pool
        'Segmentation', 'SEGNET', 'PixelClassification', ...  % semantic segmentation networks
        'reach_exact', ...                             % exact-star NNCS (memory/time blow-up)
        'comprehensive_network', ...                   % large multi-layer soundness net
        '_gpu_' ...                                    % GPU tests (no GPU on CI runners)
    };

    suite = TestSuite.fromFolder(testsDir, 'IncludingSubfolders', true);
    names = string({suite.Name});
    isHeavy = false(1, numel(names));
    for p = 1:numel(highMemPatterns)
        isHeavy = isHeavy | contains(names, highMemPatterns{p}, 'IgnoreCase', true);
    end

    switch lower(mode)
        case 'standard', sel = suite(~isHeavy);
        case 'highmem',  sel = suite(isHeavy);
        case 'all',      sel = suite;
        otherwise, error('run_shard:mode', "mode must be 'standard', 'highmem', or 'all'.");
    end

    % Deterministic round-robin partition (interleaves heavy/light tests for balance).
    sel = sel(shardIdx:numShards:end);

    fprintf('== run_shard %d/%d mode=%s: running %d of %d tests (%d heavy %s) ==\n', ...
        shardIdx, numShards, mode, numel(sel), numel(suite), nnz(isHeavy), ...
        ternary(strcmpi(mode,'standard'),'excluded','total'));

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

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
