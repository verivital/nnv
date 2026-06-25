function tests = test_nn4sys_lindex
%TEST_NN4SYS_LINDEX  Tests for the nn4sys lindex early route (sound batched-IBP + adaptive box
%   subdivision decider, nn4sys_lindex_decide.py, dispatched from run_vnncomp_instance.m). The lindex
%   spec is a huge OR over ultra-thin X-boxes; NNV's general star/LP reach times out (lindex_200 = 444s),
%   so the early route decides it via sound IBP (exact at a point -> subdividing the thin box is sound
%   AND tight). The decider is SOUND-OR-UNKNOWN with CRASH-SAFE exit codes (10=sat/11=unsat/else=unknown,
%   so a Python crash degrades to unknown, never a wrong unsat).
%
%   Groups:
%     * always-run: the decider script ships next to the runner.
%     * decide correctness (needs onnx+onnxruntime python AND the nn4sys benchmark asset): lindex_200 and
%       a lindex_deep instance -> unsat via the early route, in well under the 30s budget.
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    here = fileparts(mfilename('fullpath'));
    sub = fullfile(here, '..', '..', 'examples', 'Submission', 'VNN_COMP2026');
    tc.TestData.addedPath = '';
    if isfolder(sub) && ~contains(path, sub)
        addpath(sub); tc.TestData.addedPath = sub;
    end
    tc.TestData.sub = sub;
end

function teardownOnce(tc)
    if isfield(tc.TestData, 'addedPath') && ~isempty(tc.TestData.addedPath)
        rmpath(tc.TestData.addedPath);
    end
end

% ---------- always-run ----------

function test_decider_script_present(tc)
    verifyTrue(tc, isfile(fullfile(tc.TestData.sub, 'nn4sys_lindex_decide.py')), ...
        'nn4sys_lindex_decide.py must ship next to run_vnncomp_instance.m');
end

% ---------- decide correctness (needs python + benchmark asset) ----------

function test_lindex_unsat(tc)
    [onnx, vnn] = assumeLindexAsset(tc, 'lindex.onnx', 'lindex_200.vnnlib');
    of = [tempname '.txt']; c = onCleanup(@() delete(of)); %#ok<NASGU>
    st = run_vnncomp_instance('nn4sys', onnx, vnn, of);
    verifyEqual(tc, st, 1, 'lindex_200 must be soundly UNSAT via the early route');
    verifyTrue(tc, contains(fileread(of), 'unsat'), 'result file must say unsat');
end

function test_lindex_deep_unsat(tc)
    [onnx, vnn] = assumeLindexAsset(tc, 'lindex_deep.onnx', 'lindex_400.vnnlib');
    of = [tempname '.txt']; c = onCleanup(@() delete(of)); %#ok<NASGU>
    st = run_vnncomp_instance('nn4sys', onnx, vnn, of);
    verifyEqual(tc, st, 1, 'lindex_deep/lindex_400 must be soundly UNSAT (deeper net, IBP looser -> more subdivision)');
end

% ---------- helpers ----------

function [onnx, vnn] = assumeLindexAsset(tc, onnxName, vnnName)
    assumeFalse(tc, isempty(resolve_ort_python()), 'no python imports onnx+onnxruntime -- skip');
    here = fileparts(mfilename('fullpath'));
    base = fullfile(here, '..', '..', '..', '..', '..', 'vnncomp2026_benchmarks', 'benchmarks', 'nn4sys', '2.0');
    onnx = fullfile(base, 'onnx', onnxName);
    vnn  = fullfile(base, 'vnnlib', vnnName);
    if ~isfile(onnx) && isfile([onnx '.gz']), gunzip([onnx '.gz']); end
    if ~isfile(vnn)  && isfile([vnn '.gz']),  gunzip([vnn '.gz']);  end
    assumeTrue(tc, isfile(onnx) && isfile(vnn), 'nn4sys lindex benchmark asset not present -- skip');
end

function py = resolve_ort_python()
    cands = {strtrim(getenv('NNV_ORT_PYTHON')), fullfile(getenv('HOME'), 'taylor_venv', 'bin', 'python'), 'python3', 'python'};
    py = '';
    for i = 1:numel(cands)
        c = cands{i}; if isempty(c), continue; end
        [st, ~] = system(sprintf('"%s" -c "import onnx, onnxruntime, numpy"', c));
        if st == 0, py = c; return; end
    end
end
