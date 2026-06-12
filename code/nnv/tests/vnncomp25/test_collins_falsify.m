function tests = test_collins_falsify
%TEST_COLLINS_FALSIFY  Tests for the collins_aerospace_benchmark falsification path.
%   collins_aerospace_benchmark is handled ENTIRELY by collins_falsify.py (sound
%   sat-or-unknown falsification through onnxruntime; UNSAT is out of scope -- the
%   YOLOv5 Detect head defeats MATLAB's importer and set-based reach at 640x640x3 is
%   infeasible). run_vnncomp_instance routes the category to the script BEFORE any
%   importNetworkFromONNX call, so the old 'Unsupported Class of Layer' failure (and
%   the older invalid-SAT CHW/HWC mix-up, a guaranteed -150) cannot recur.
%
%   Test groups:
%     * routing contract (always runs): the dispatcher must route collins to the
%       falsifier and complete with 'unknown' on bogus inputs -- no importer error.
%     * falsifier correctness (gated, like test_validate_witness_onnx, on
%       python+onnxruntime AND the local vnncomp2026 benchmark files):
%         - center-consistency gate passes on a real spec (proves the HWC mapping);
%         - the known-falsifiable delta=0.1 instance reproduces SAT end-to-end
%           through run_vnncomp_instance, with a well-formed counterexample.
%
%   The gated tests run python on ~140MB (decompressed) specs: the center check
%   takes ~10-30s; the end-to-end SAT reproduction took ~4 min on the dev box
%   (2984 onnxruntime forwards). Set NNV_COLLINS_BUDGET (seconds) to bound the
%   SAT search (the SAT test sets 600).
%
% Run: results = runtests('test_collins_falsify')
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    here = fileparts(mfilename('fullpath'));
    sub = fullfile(here, '..', '..', 'examples', 'Submission', 'VNN_COMP2025');
    tc.TestData.subdir = sub;
    tc.TestData.addedPath = false;
    if isfolder(sub) && ~contains(path, sub)
        addpath(sub);
        tc.TestData.addedPath = true;
    end
    tc.assumeTrue(exist('run_vnncomp_instance', 'file') == 2, ...
        'run_vnncomp_instance not on path; add examples/Submission/VNN_COMP2025');
    tc.TestData.script = fullfile(sub, 'collins_falsify.py');
    % Local benchmark files (kept .gz on disk; the script gunzips transparently).
    % <repo>/code/nnv/tests/vnncomp25 -> workspace root is 5 levels up.
    bench = fullfile(here, '..', '..', '..', '..', '..', ...
        'vnncomp2026_benchmarks', 'benchmarks', 'collins_aerospace_benchmark', '1.0');
    tc.TestData.onnx   = fullfile(bench, 'onnx', 'yolov5nano_LRelu_640.onnx.gz');
    tc.TestData.spec01  = fullfile(bench, 'vnnlib', 'img_12761_perturbed_bbox_0_delta_0.1.vnnlib.gz');
    tc.TestData.spec001 = fullfile(bench, 'vnnlib', 'img_14421_perturbed_bbox_3_delta_0.001.vnnlib.gz');
end

function teardownOnce(tc)
    if isfield(tc.TestData, 'addedPath') && tc.TestData.addedPath
        rmpath(tc.TestData.subdir);
    end
end

% --------------------------------------------------------------------------
% routing contract (environment-independent)
% --------------------------------------------------------------------------

function test_falsifier_script_exists(tc)
    verifyTrue(tc, isfile(tc.TestData.script), ...
        'collins_falsify.py must ship next to run_vnncomp_instance.m');
end

function test_dispatch_routes_to_falsifier_not_importer(tc)
    % With bogus inputs the collins branch must COMPLETE with 'unknown' -- the
    % python call fails to parse a non-existent vnnlib and the wrapper degrades
    % safely. Before this path existed, dispatch went to importNetworkFromONNX and
    % died with 'Unsupported Class of Layer' on the Detect-head custom layers; any
    % error here (importer or otherwise) is a routing regression.
    onnx = [tempname '.onnx'];        % does not exist
    vnnlib = [tempname '.vnnlib'];    % does not exist
    outf = [tempname '.txt'];
    cu = onCleanup(@() delete_if(outf));
    try
        status = run_vnncomp_instance('collins_aerospace_benchmark', onnx, vnnlib, outf);
    catch ME
        verifyFail(tc, sprintf(['collins dispatch must not error (was: importer ' ...
            '''Unsupported Class of Layer''); got %s: %s'], ME.identifier, ME.message));
        return;
    end
    verifyEqual(tc, status, 2, 'bogus inputs must yield unknown (2), never sat/unsat');
    txt = strtrim(fileread(outf));
    verifyEqual(tc, txt, 'unknown', 'output file must say unknown');
end

% --------------------------------------------------------------------------
% falsifier correctness (needs python+onnxruntime and the local benchmarks)
% --------------------------------------------------------------------------

function test_center_consistency_proves_hwc_mapping(tc)
    % --check-center-only: parse the real delta=0.001 spec, run the box center
    % through onnxruntime with the HWC-flat -> CHW mapping, and confirm the output
    % property does NOT hold there (objectness inside its band). Under a wrong
    % (CHW-read) mapping the center looks violated and the script refuses; the
    % CENTER_CONSISTENT marker is therefore positive proof of the mapping.
    assumeFalsifierEnv(tc);
    cmd = sprintf('%s "%s" "%s" "%s" "%s" 60 --check-center-only', resolve_python(), ...
        tc.TestData.script, tc.TestData.onnx, tc.TestData.spec001, [tempname '.csv']);
    [st, out] = system(cmd);
    verifyTrue(tc, contains(out, 'CENTER_CONSISTENT'), ...
        sprintf('center-consistency must pass on the real spec; output:\n%s', out));
    verifyEqual(tc, st, 2, '--check-center-only is diagnostic: must exit 2 (unknown), not claim sat');
end

function test_delta01_sat_reproduces_end_to_end(tc)
    % The delta=0.1 instance is known falsifiable (objectness 0.8809 -> below the
    % 0.792 band edge in one FD-gradient step). Run the FULL dispatcher and expect
    % a sat verdict plus a well-formed counterexample (X_0..X_1228799, Y_0..).
    assumeFalsifierEnv(tc);
    setenv('NNV_COLLINS_BUDGET', '600');
    cu0 = onCleanup(@() setenv('NNV_COLLINS_BUDGET', ''));
    outf = [tempname '.txt'];
    cu = onCleanup(@() delete_if(outf));
    status = run_vnncomp_instance('collins_aerospace_benchmark', ...
        tc.TestData.onnx, tc.TestData.spec01, outf);
    verifyEqual(tc, status, 0, 'delta=0.1 must falsify (sat)');
    txt = fileread(outf);
    verifyTrue(tc, startsWith(strtrim(txt), 'sat'), 'output file must start with sat');
    verifyTrue(tc, contains(txt, '(X_0 '), 'counterexample must assign X_0');
    verifyTrue(tc, contains(txt, '(X_1228799 '), 'counterexample must assign all 1228800 inputs');
    verifyTrue(tc, contains(txt, '(Y_0 '), 'counterexample must assign outputs');
end

% --------------------------------------------------------------------------
% helpers
% --------------------------------------------------------------------------

function assumeFalsifierEnv(tc)
    % Same python resolution as the wrapper (validate_witness_onnx idiom), and the
    % local 2026 benchmark files must exist. Skip -- never fail -- otherwise.
    [st, ~] = system([resolve_python() ' -c "import onnxruntime, numpy"']);
    assumeEqual(tc, st, 0, 'python/onnxruntime unavailable -- skipping falsifier tests');
    assumeTrue(tc, isfile(tc.TestData.onnx) && isfile(tc.TestData.spec01) ...
        && isfile(tc.TestData.spec001), ...
        'vnncomp2026 collins benchmark files not found locally -- skipping');
end

function py = resolve_python()
    py = 'python';
    try
        pe = pyenv;
        if ~isempty(pe.Executable) && isfile(pe.Executable)
            py = ['"' char(pe.Executable) '"'];
        end
    catch
    end
end

function delete_if(f)
    if exist(f, 'file'), delete(f); end
end
