function run_shard(shardIdx, numShards, mode)
% RUN_SHARD  Run a deterministic, duration-balanced 1-of-N shard of the NNV test suite.
%
%   run_shard(k, N)             % standard shard k of N
%   run_shard(k, N, 'standard') % same
%   run_shard(k, N, 'highmem')  % ONLY the anticipated-heavy tests (big-RAM/self-hosted/ACCRE)
%   run_shard(k, N, 'all')      % everything CI can run (excludes only GPU tests)
%
% For GitHub Actions strategy.matrix via matlab-actions/run-command:
%   cd('code/nnv'); install; addpath(genpath(fullfile(pwd,'tests')));
%   run_shard(${{ matrix.shard }}, ${{ N }}, 'standard');
%
% - Sharding is DURATION-BALANCED (LPT bin-packing over test_durations.csv) to minimize
%   wall-clock; falls back to round-robin if the CSV is missing.
% - Writes JUnit XML + a ProgressPlugin START/DONE trace to code/nnv/test-results-ci/
%   (uploaded as a CI artifact). It does NOT assertSuccess: failing tests are recorded and
%   gated by the report job (ci_report.py vs ci_allowed_failures.txt). A shard fails only if
%   MATLAB itself crashes -- the progress trace then names the offending test.

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

    % Excluded from CI shards (cannot run on a GitHub-hosted runner):
    %   GPU tests -- no GPU on the runners (error/incomplete there).
    gpuPatterns = {'_gpu_', 'toGPU'};
    %   CONFIRMED to crash MATLAB on the ~16 GB runner (exit 255 / OOM). Excluded so the
    %   shard's other results survive; run these via the high-mem fallback (highmem mode,
    %   self-hosted/ACCRE) or locally. Populate (exact full Names) from the report's crash list.
    confirmedCrashPatterns = { };

    % Tests that merely FAIL (incl. WIP transformer/SLM, tutorial env-diffs) are NOT excluded:
    % they run, results land in JUnit, and the report job classifies them vs
    % ci_allowed_failures.txt (known = visible + non-blocking; new = red).

    % Anticipated-heavy (memory/time) tests. EMPIRICALLY these flake on the 16 GB runner:
    % including them caused NONDETERMINISTIC crashes (exit 255) AND hangs (a shard ran 30+
    % min while peers finished in ~10). So they are EXCLUDED from the standard shards and run
    % via mode='highmem' on a big-RAM/self-hosted runner or ACCRE (or locally). To recover
    % specific ones into standard later, pinpoint the few flaky tests and narrow this list.
    anticipatedHeavyPatterns = {'VolumeStar', 'Conv3D', 'AveragePooling3D', ...
        'Segmentation', 'SEGNET', 'PixelClassification', 'reach_exact', 'comprehensive_network'};

    ciExclude = [gpuPatterns, anticipatedHeavyPatterns, confirmedCrashPatterns];

    suite = TestSuite.fromFolder(testsDir, 'IncludingSubfolders', true);
    names = string({suite.Name});

    switch lower(mode)
        case 'standard', sel = suite(~matchAny(names, ciExclude));
        case 'highmem',  sel = suite(matchAny(names, anticipatedHeavyPatterns) & ~matchAny(names, gpuPatterns));
        case 'all',      sel = suite(~matchAny(names, gpuPatterns));
        otherwise, error('run_shard:mode', "mode must be 'standard', 'highmem', or 'all'.");
    end

    fprintf('== run_shard %d/%d mode=%s: %d of %d total (excluded: %d gpu, %d heavy, %d confirmed-crash) ==\n', ...
        shardIdx, numShards, mode, numel(sel), numel(suite), ...
        nnz(matchAny(names, gpuPatterns)), nnz(matchAny(names, anticipatedHeavyPatterns)), nnz(matchAny(names, confirmedCrashPatterns)));

    % Duration-balanced (LPT bin-packing) sharding: assign each test heaviest-first to the
    % currently-lightest shard. Deterministic, so every shard computes the same global
    % assignment and selects its own bin. Round-robin fallback if test_durations.csv is
    % missing. (Regenerate the CSV from a full-suite run:
    %   writetable(sortrows(T(:,{'Name','Duration'}),'Duration','descend'), test_durations.csv).)
    selNames = string({sel.Name});
    durFile = fullfile(testsDir, 'test_durations.csv');
    if isfile(durFile) && numShards > 1
        DT = readtable(durFile, 'TextType', 'string');
        dmap = containers.Map('KeyType','char','ValueType','double');
        for i = 1:height(DT), dmap(char(DT.Name(i))) = DT.Duration(i); end
        defDur = 1.0;   % default seconds for new/unknown tests
        durs = zeros(1, numel(selNames));
        for i = 1:numel(selNames)
            k = char(selNames(i));
            if isKey(dmap, k), durs(i) = dmap(k); else, durs(i) = defDur; end
        end
        [~, order] = sort(durs, 'descend');   % stable on ties -> deterministic
        binTot = zeros(1, numShards); binOf = zeros(1, numel(durs));
        for j = order
            [~, b] = min(binTot);
            binOf(j) = b; binTot(b) = binTot(b) + durs(j);
        end
        sel = sel(binOf == shardIdx);
        fprintf('  (duration-balanced: this shard est %.0fs of %.0fs total / %d shards)\n', ...
            binTot(shardIdx), sum(durs), numShards);
    else
        sel = sel(shardIdx:numShards:end);   % round-robin fallback
    end

    outDir = fullfile(fileparts(testsDir), 'test-results-ci');   % code/nnv/test-results-ci
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    xmlFile = fullfile(outDir, sprintf('results-%s-shard%02d-of%02d.xml', lower(mode), shardIdx, numShards));
    fprintf('JUnit report -> %s\n', xmlFile);

    if isempty(sel)
        fprintf('No tests selected for this shard; writing empty report and exiting success.\n');
        fid = fopen(xmlFile, 'w'); fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n<testsuites/>\n'); fclose(fid);
        return;
    end

    progressFile = fullfile(outDir, sprintf('progress-%s-shard%02d-of%02d.txt', lower(mode), shardIdx, numShards));
    runner = TestRunner.withTextOutput('OutputDetail', Verbosity.Concise);
    runner.addPlugin(XMLPlugin.producingJUnitFormat(xmlFile));
    runner.addPlugin(ProgressPlugin(progressFile));   % self-identifying crash log (last START w/o DONE = crasher)
    results = runner.run(sel);

    fprintf('== shard %d/%d done: %d passed, %d failed, %d incomplete (failures gated by report job) ==\n', ...
        shardIdx, numShards, nnz([results.Passed]), nnz([results.Failed]), nnz([results.Incomplete]));
end

function tf = matchAny(names, patterns)
    tf = false(1, numel(names));
    for p = 1:numel(patterns)
        tf = tf | contains(names, patterns{p}, 'IgnoreCase', true);
    end
end
