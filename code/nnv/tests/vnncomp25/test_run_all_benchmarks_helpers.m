function tests = test_run_all_benchmarks_helpers
% test_run_all_benchmarks_helpers
% Function-based unit tests for the PURE helpers of run_all_benchmarks.m that
% need no network / parpool: instances.csv parsing + status-code stringification.
%
% NOTE ON TESTABILITY (see report): the production logic lives in LOCAL functions
% `pick_instances` / `ensure_decompressed` / `status_to_str` INSIDE
% run_all_benchmarks.m, which are not reachable from a test. This test exercises
% the behavior-preserving extracted copies:
%     examples/Submission/VNN_COMP2026/vnncomp_pick_instances.m
%     examples/Submission/VNN_COMP2026/vnncomp_status_to_str.m
% Recommended follow-up refactor: have run_all_benchmarks.m delegate to these so
% the tested code IS the production code.
%
% Run: results = runtests('test_run_all_benchmarks_helpers')
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Make the extracted helpers visible on the path.
    here = fileparts(mfilename('fullpath'));
    sub = fullfile(here, '..', '..', 'examples', 'Submission', 'VNN_COMP2026');
    testCase.TestData.subdir = sub;
    testCase.TestData.addedPath = false;
    if isfolder(sub) && ~contains(path, sub)
        addpath(sub);
        testCase.TestData.addedPath = true;
    end
end

function teardownOnce(testCase)
    if isfield(testCase.TestData, 'addedPath') && testCase.TestData.addedPath
        rmpath(testCase.TestData.subdir);
    end
end

% --------------------------------------------------------------------------
% status_to_str 0/1/2 + negatives + fallthrough
% --------------------------------------------------------------------------
function test_status_to_str_mapping(testCase)
    verifyEqual(testCase, vnncomp_status_to_str(0),  'sat');
    verifyEqual(testCase, vnncomp_status_to_str(1),  'unsat');
    verifyEqual(testCase, vnncomp_status_to_str(2),  'unknown');
    verifyEqual(testCase, vnncomp_status_to_str(-1), 'error');
    verifyEqual(testCase, vnncomp_status_to_str(-2), 'missing');
    verifyEqual(testCase, vnncomp_status_to_str(-3), 'decompress_failed');
    % unrecognized code -> code_<n>
    verifyEqual(testCase, vnncomp_status_to_str(7),   'code_7');
    verifyEqual(testCase, vnncomp_status_to_str(-99), 'code_-99');
end

% --------------------------------------------------------------------------
% pick_instances 'first' returns exactly ONE resolvable pair, './'-stripped
% --------------------------------------------------------------------------
function test_pick_instances_first(testCase)
    [root, sub] = make_bench(testCase, { ...
        './onnx/a.onnx', './vnnlib/p1.vnnlib'; ...
        './onnx/b.onnx', './vnnlib/p2.vnnlib'; ...
        './onnx/c.onnx', './vnnlib/p3.vnnlib'}, ...
        {'onnx/a.onnx','vnnlib/p1.vnnlib','onnx/b.onnx','vnnlib/p2.vnnlib', ...
         'onnx/c.onnx','vnnlib/p3.vnnlib'});
    inst = fullfile(root, sub, 'instances.csv');
    pairs = vnncomp_pick_instances(inst, root, sub, 'first');

    verifyEqual(testCase, numel(pairs), 1, 'first -> exactly one instance');
    verifyEqual(testCase, pairs{1}.onnx_rel,   'onnx/a.onnx', './ stripped from onnx_rel');
    verifyEqual(testCase, pairs{1}.vnnlib_rel, 'vnnlib/p1.vnnlib', './ stripped from vnnlib_rel');
    verifyTrue(testCase, isfile(pairs{1}.onnx),   'resolved onnx path exists');
    verifyTrue(testCase, isfile(pairs{1}.vnnlib), 'resolved vnnlib path exists');
end

% --------------------------------------------------------------------------
% pick_instances 'all' returns EVERY resolvable pair, in CSV order
% --------------------------------------------------------------------------
function test_pick_instances_all(testCase)
    [root, sub] = make_bench(testCase, { ...
        'onnx/a.onnx', 'vnnlib/p1.vnnlib'; ...
        'onnx/b.onnx', 'vnnlib/p2.vnnlib'; ...
        'onnx/c.onnx', 'vnnlib/p3.vnnlib'}, ...
        {'onnx/a.onnx','vnnlib/p1.vnnlib','onnx/b.onnx','vnnlib/p2.vnnlib', ...
         'onnx/c.onnx','vnnlib/p3.vnnlib'});
    inst = fullfile(root, sub, 'instances.csv');
    pairs = vnncomp_pick_instances(inst, root, sub, 'all');

    verifyEqual(testCase, numel(pairs), 3, 'all -> every resolvable instance');
    verifyEqual(testCase, pairs{1}.onnx_rel, 'onnx/a.onnx');
    verifyEqual(testCase, pairs{2}.onnx_rel, 'onnx/b.onnx');
    verifyEqual(testCase, pairs{3}.onnx_rel, 'onnx/c.onnx');
    verifyEqual(testCase, pairs{2}.vnnlib_rel, 'vnnlib/p2.vnnlib');
end

% --------------------------------------------------------------------------
% pick_instances skips rows whose files do NOT exist (unresolvable)
% --------------------------------------------------------------------------
function test_pick_instances_skips_unresolvable(testCase)
    % Only create files for rows 1 and 3; row 2 references missing files.
    [root, sub] = make_bench(testCase, { ...
        'onnx/a.onnx', 'vnnlib/p1.vnnlib'; ...
        'onnx/missing.onnx', 'vnnlib/missing.vnnlib'; ...
        'onnx/c.onnx', 'vnnlib/p3.vnnlib'}, ...
        {'onnx/a.onnx','vnnlib/p1.vnnlib','onnx/c.onnx','vnnlib/p3.vnnlib'});
    inst = fullfile(root, sub, 'instances.csv');

    % 'all' returns only the two resolvable rows
    pairs = vnncomp_pick_instances(inst, root, sub, 'all');
    verifyEqual(testCase, numel(pairs), 2, 'unresolvable middle row skipped');
    verifyEqual(testCase, pairs{1}.onnx_rel, 'onnx/a.onnx');
    verifyEqual(testCase, pairs{2}.onnx_rel, 'onnx/c.onnx');

    % 'first' still returns the FIRST resolvable row (row 1 here)
    pairs1 = vnncomp_pick_instances(inst, root, sub, 'first');
    verifyEqual(testCase, numel(pairs1), 1);
    verifyEqual(testCase, pairs1{1}.onnx_rel, 'onnx/a.onnx');
end

% --------------------------------------------------------------------------
% pick_instances returns {} when nothing is resolvable
% --------------------------------------------------------------------------
function test_pick_instances_none_resolvable(testCase)
    [root, sub] = make_bench(testCase, { ...
        'onnx/x.onnx', 'vnnlib/x.vnnlib'}, ...
        {});   % create NO files
    inst = fullfile(root, sub, 'instances.csv');
    pairs = vnncomp_pick_instances(inst, root, sub, 'all');
    verifyTrue(testCase, isempty(pairs), 'no resolvable instances -> empty');
end

% --------------------------------------------------------------------------
% Helpers
% --------------------------------------------------------------------------
function [root, sub] = make_bench(testCase, rows, files_to_create)
    % rows: N x 2 cellstr of {onnx_rel, vnnlib_rel} (timeout col auto-added).
    % files_to_create: cellstr of relative paths (under root/sub) to materialize.
    root = tempname; mkdir(root);
    sub = 'benchX';
    mkdir(fullfile(root, sub));
    mkdir(fullfile(root, sub, 'onnx'));
    mkdir(fullfile(root, sub, 'vnnlib'));
    % write instances.csv (3 columns: onnx_rel, vnnlib_rel, timeout)
    csvfile = fullfile(root, sub, 'instances.csv');
    fid = fopen(csvfile, 'w');
    for r = 1:size(rows, 1)
        fprintf(fid, '%s,%s,%d\n', rows{r,1}, rows{r,2}, 120);
    end
    fclose(fid);
    % materialize the requested files
    for k = 1:numel(files_to_create)
        f = fullfile(root, sub, files_to_create{k});
        fd = fopen(f, 'w'); fprintf(fd, 'x'); fclose(fd);
    end
    % cleanup at end of each test
    testCase.addTeardown(@() rmdir(root, 's'));
end
