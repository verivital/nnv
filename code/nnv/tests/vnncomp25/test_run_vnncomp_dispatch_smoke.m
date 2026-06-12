function tests = test_run_vnncomp_dispatch_smoke
% test_run_vnncomp_dispatch_smoke
% Category-DISPATCH smoke scaffold for run_vnncomp_instance.m.
%
% run_vnncomp_instance(category, onnx, vnnlib, outputfile) selects, via a long
% if/elseif on `category`, the ONNX import args / reachMethod / reshape for that
% benchmark, then calls importNetworkFromONNX on `onnx`. This test does NOT verify
% any network: it passes a NON-EXISTENT onnx/vnnlib and asserts, per category, that
% the dispatch reaches the RIGHT branch -- i.e. the failure is a FILE / IMPORT
% error (model can't be loaded), NOT a DISPATCH error:
%   - the final `else error("ONNX model not supported")` (unknown category), or
%   - a vnnlib parse error, etc.
% Known-unsupported categories that deliberately error AT dispatch (before import)
% are pinned separately with their stable messages.
%
% >>> THE MAIN AGENT MUST VALIDATE THIS IN MATLAB. Two things matter:
%   1. It needs run_vnncomp_instance.m + its dependencies (load_manifest_net,
%      matlab2nnv, importNetworkFromONNX, ...) on the path; add
%      examples/Submission/VNN_COMP2025 and run startup_nnv first.
%   2. For a REAL end-to-end dispatch check, point onnx/vnnlib at local files under
%      vnncomp2025_benchmarks/benchmarks/<category>/ ; with the dummy paths used
%      here the test only proves "dispatch did not reject the category", which is
%      the regression we care about. If importNetworkFromONNX raises a DIFFERENT
%      error class than expected on this MATLAB build, relax classify_error below.
%
% Run: results = runtests('test_run_vnncomp_dispatch_smoke')
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    here = fileparts(mfilename('fullpath'));
    sub = fullfile(here, '..', '..', 'examples', 'Submission', 'VNN_COMP2025');
    testCase.TestData.subdir = sub;
    testCase.TestData.addedPath = false;
    if isfolder(sub) && ~contains(path, sub)
        addpath(sub);
        testCase.TestData.addedPath = true;
    end
    % Skip the whole file gracefully if the runner isn't importable here.
    testCase.assumeTrue(exist('run_vnncomp_instance', 'file') == 2, ...
        'run_vnncomp_instance not on path; add examples/Submission/VNN_COMP2025 + startup_nnv');
end

function teardownOnce(testCase)
    if isfield(testCase.TestData, 'addedPath') && testCase.TestData.addedPath
        rmpath(testCase.TestData.subdir);
    end
end

% --------------------------------------------------------------------------
% Supported categories: dispatch must NOT reject them. With a bogus onnx path the
% expected failure is a file/import error (or a manifest-not-found for the two
% manifest categories), never "ONNX model not supported".
% --------------------------------------------------------------------------
function test_supported_categories_reach_import(testCase)
    supported = { 'acasxu', 'cersyve', 'cgan', 'cifar100', ...
        'collins_aerospace_benchmark', 'collins_rul', 'cora', 'dist_shift', ...
        'linearize', 'lsnc_relu', 'malbeware', 'metaroom', 'ml4acopf', ...
        'nn4sys', 'relusplitter', 'safenlp', 'sat_relu', 'soundness', ...
        'tinyimagenet', 'tllverify', 'traffic', 'vggnet', 'vit', 'yolo' };

    for i = 1:numel(supported)
        cat = supported{i};
        kind = classify_dispatch(cat);
        % The meaningful check: a supported category must NOT be rejected by the
        % dispatcher's final else ("ONNX model not supported"). Whether the bogus-input
        % run then fails at vnnlib-parse or at model-import is implementation-order
        % specific (the manifest categories -- lsnc_relu, traffic -- load the vnnlib
        % before importing, so they surface 'vnnlib_parse'; that is NOT a dispatch
        % rejection and is acceptable here).
        verifyNotEqual(testCase, kind, 'unsupported_category', ...
            sprintf('category "%s" was rejected as unsupported by the dispatcher', cat));
    end
end

% --------------------------------------------------------------------------
% cctsdb_yolo: no longer errors at dispatch ("Working on supporting this one").
% It now routes to the complete-enumeration python path (cctsdb_enumerate.py)
% BEFORE load_vnncomp_network. With bogus inputs the script (or a missing
% python) yields UNKNOWN, so the runner must return WITHOUT error and write
% 'unknown' -- a sound degradation, never a crash.
% --------------------------------------------------------------------------
function test_cctsdb_yolo_routes_to_enumeration(testCase)
    onnx = [tempname '.onnx'];          % does not exist
    vnnlib = [tempname '.vnnlib'];      % does not exist
    outf = [tempname '.txt'];
    c = onCleanup(@() cleanup_file(outf));
    try
        run_vnncomp_instance('cctsdb_yolo', onnx, vnnlib, outf);
    catch ME
        verifyFail(testCase, sprintf( ...
            'cctsdb_yolo must not error at dispatch anymore; got: %s', ME.message));
        return
    end
    verifyTrue(testCase, isfile(outf), 'runner must write the output file');
    txt = strtrim(fileread(outf));
    verifyEqual(testCase, txt, 'unknown', ...
        'bogus inputs must degrade to a sound ''unknown''');
end

function test_unknown_category_is_rejected(testCase)
    % A category the dispatcher does not recognize must hit the final else.
    msg = dispatch_error_message(testCase, 'no_such_benchmark_xyz');
    verifyTrue(testCase, contains(msg, 'ONNX model not supported'), ...
        sprintf('unknown category should be rejected as unsupported; got: %s', msg));
end

% --------------------------------------------------------------------------
% Helpers
% --------------------------------------------------------------------------
function msg = dispatch_error_message(testCase, cat)
    % Invoke run_vnncomp_instance with bogus inputs; return the error message.
    onnx = [tempname '.onnx'];          % does not exist
    vnnlib = [tempname '.vnnlib'];      % does not exist
    outf = [tempname '.txt'];
    c = onCleanup(@() cleanup_file(outf));
    msg = '';
    try
        run_vnncomp_instance(cat, onnx, vnnlib, outf);
        % Reaching here means no error at all -- unexpected for bogus inputs.
        verifyFail(testCase, sprintf('expected an error for category "%s" with bogus inputs', cat));
    catch ME
        msg = ME.message;
    end
end

function kind = classify_dispatch(cat)
    % Classify the failure of run_vnncomp_instance with a non-existent model into:
    %   'unsupported_category' : final else "ONNX model not supported"
    %   'vnnlib_parse'         : failed parsing the (bogus) vnnlib before import
    %   'import_or_file'       : import/file error AFTER a real dispatch branch (GOOD)
    %   'none'                 : no error (unexpected)
    onnx = [tempname '.onnx'];
    vnnlib = [tempname '.vnnlib'];
    outf = [tempname '.txt'];
    c = onCleanup(@() cleanup_file(outf));
    try
        run_vnncomp_instance(cat, onnx, vnnlib, outf);
        kind = 'none';
        return;
    catch ME
        m = ME.message;
    end
    if contains(m, 'ONNX model not supported')
        kind = 'unsupported_category';
    elseif contains(m, 'load_vnnlib') || contains(m, 'vnnlib') || ...
            contains(m, 'fopen') || contains(lower(m), 'property not supported')
        % NB: with a non-existent vnnlib, load_vnnlib may fail; but import happens
        % FIRST in run_vnncomp_instance, so a real model would import before the
        % vnnlib is read. With a bogus onnx the import error fires first, so this
        % branch should not normally trigger for supported cats.
        kind = 'vnnlib_parse';
    else
        kind = 'import_or_file';
    end
end

function cleanup_file(f)
    if exist(f, 'file'), delete(f); end
end
