function tests = test_authoritative_witness_gate
%TEST_AUTHORITATIVE_WITNESS_GATE  Tests for authoritative_witness_gate.m -- the run_all_benchmarks
%   (dev-sweep) twin of execute.py's authoritative SAT-witness gate. The helper re-validates an
%   emitted `sat` result file against the AUTHORITATIVE vnnlib parser + a real-onnx onnxruntime
%   forward (validate_witness_authoritative.py) and returns 'valid' / 'spurious' / 'cant'. Wiring it
%   into run_all_benchmarks is what makes a sweep scorecard PER-WITNESS trustworthy (vs only the
%   consensus-level check the 2026-06-24 sweep had).
%
%   Two groups:
%     * SOUNDNESS / fail-open contract (ALWAYS runs, no deps): every degenerate input -- missing
%       result file, non-sat verdict, unparseable onnx/vnnlib -- maps to 'cant'. The gate NEVER
%       returns 'valid'/'spurious' when it cannot actually check, and the verdict is always one of
%       the three known strings. This is the property that keeps the gate sound-or-keep.
%     * replay correctness (runs only if a python imports onnx+onnxruntime+vnnlib AND the acasxu
%       benchmark asset is present): a made-up witness on a real spec -> 'spurious'; a STACKED
%       two-network witness -> 'cant' (the KEYSTONE -- the gate provably CANNOT catch the iso/mono
%       case, which is exactly why run_vnncomp_instance sound-downgrades the multinet branch at the
%       source; without that downgrade those 50 spurious monotonic_acasxu sats would be -7500).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    % authoritative_witness_gate lives in examples/Submission/VNN_COMP2026. CI shards tests across
    % matrix jobs and does not guarantee that dir is on the path, so add it here to be self-sufficient.
    here = fileparts(mfilename('fullpath'));
    sub = fullfile(here, '..', '..', 'examples', 'Submission', 'VNN_COMP2026');
    tc.TestData.addedPath = '';
    if isfolder(sub) && ~contains(path, sub)
        addpath(sub); tc.TestData.addedPath = sub;
    end
end

function teardownOnce(tc)
    if isfield(tc.TestData, 'addedPath') && ~isempty(tc.TestData.addedPath)
        rmpath(tc.TestData.addedPath);
    end
end

% ---------- soundness / fail-open contract (environment-independent) ----------

function test_missing_result_file_is_cant(tc)
    v = authoritative_witness_gate('no.onnx', 'no.vnnlib', fullfile(tempdir, 'nope_xyz_123.txt'));
    verifyEqual(tc, v, 'cant', 'a missing result file must fail-open (keep the verdict)');
end

function test_nonsat_file_is_cant(tc)
    f = [tempname '.txt']; fid = fopen(f, 'w'); fprintf(fid, 'unsat\n'); fclose(fid);
    c = onCleanup(@() delete(f)); %#ok<NASGU>
    v = authoritative_witness_gate('no.onnx', 'no.vnnlib', f);
    verifyEqual(tc, v, 'cant', 'a non-sat result must never be judged spurious/valid');
end

function test_bad_paths_with_sat_file_is_cant(tc)
    % A syntactically real sat witness but bogus onnx/vnnlib -> the script cannot parse/forward ->
    % the gate must fail-open, NEVER return spurious/valid on an uncheckable input.
    f = [tempname '.txt']; fid = fopen(f, 'w'); fprintf(fid, 'sat\n(X_0 0.5)\n(X_1 0.5)\n'); fclose(fid);
    c = onCleanup(@() delete(f)); %#ok<NASGU>
    v = authoritative_witness_gate(fullfile(tempdir, 'no_model_xyz.onnx'), ...
                                   fullfile(tempdir, 'no_spec_xyz.vnnlib'), f);
    verifyEqual(tc, v, 'cant', 'unparseable onnx/vnnlib must fail-open');
end

function test_verdict_is_always_known_string(tc)
    f = [tempname '.txt']; fid = fopen(f, 'w'); fprintf(fid, 'sat\n(X_0 0.5)\n'); fclose(fid);
    c = onCleanup(@() delete(f)); %#ok<NASGU>
    v = authoritative_witness_gate('x.onnx', 'y.vnnlib', f);
    verifyTrue(tc, ismember(v, {'valid', 'spurious', 'cant'}), 'verdict must be one of the three');
end

% ---------- replay correctness (needs python+onnx+onnxruntime+vnnlib AND the acasxu asset) ----------

function test_madeup_witness_spurious(tc)
    [onnx, vnnlib] = assumeAcasxuAsset(tc);
    f = writeSat([0.5 0.5 0.5 0.5 0.5]);   % a made-up 5-vector that does not satisfy the real spec
    c = onCleanup(@() delete(f)); %#ok<NASGU>
    v = authoritative_witness_gate(onnx, vnnlib, f);
    verifyEqual(tc, v, 'spurious', 'a made-up witness must be caught -> sound sat->unknown downgrade');
end

function test_stacked_multinet_witness_is_cant(tc)
    % KEYSTONE: a STACKED two-network witness (10 values for a 5-input acasxu net) cannot be replayed
    % by the single-InferenceSession authoritative gate -> 'cant'. This proves the gate will NOT catch
    % the iso/mono (monotonic_acasxu) spurious sats, which is precisely why run_vnncomp_instance's
    % multinet branch sound-downgrades to unknown at the SOURCE (preventing the would-be -7500).
    [onnx, vnnlib] = assumeAcasxuAsset(tc);
    f = writeSat([0.1 0.1 0.1 0.1 0.1 0.2 0.2 0.2 0.2 0.2]);   % 10 values vs a 5-input net
    c = onCleanup(@() delete(f)); %#ok<NASGU>
    v = authoritative_witness_gate(onnx, vnnlib, f);
    verifyEqual(tc, v, 'cant', 'a stacked multinet witness is un-replayable -> gate fail-opens');
end

% ---------- helpers ----------

function [onnx, vnnlib] = assumeAcasxuAsset(tc)
    % Skip (not fail) when the witness-gate python or the acasxu benchmark asset is unavailable, so a
    % CI box without onnx/ort/vnnlib OR without the benchmarks repo stays green.
    assumeFalse(tc, isempty(resolve_gate_python()), 'no python imports onnx+onnxruntime+vnnlib -- skip');
    here = fileparts(mfilename('fullpath'));
    base = fullfile(here, '..', '..', '..', '..', '..', 'vnncomp2026_benchmarks', 'benchmarks', 'acasxu_2023');
    onnx = '';
    for v = {'2.0', '1.0'}
        cand = fullfile(base, v{1}, 'onnx', 'ACASXU_run2a_1_1_batch_2000.onnx');
        if isfile(cand) || isfile([cand '.gz']), onnx = cand; vdir = fullfile(base, v{1}, 'vnnlib'); break; end
    end
    assumeFalse(tc, isempty(onnx), 'acasxu benchmark asset not present -- skip replay test');
    vnnlib = fullfile(vdir, 'prop_3.vnnlib');
    assumeTrue(tc, isfile(vnnlib) || isfile([vnnlib '.gz']), 'acasxu prop_3.vnnlib not present -- skip');
end

function py = resolve_gate_python()
    cands = {strtrim(getenv('NNV_ORT_PYTHON')), fullfile(getenv('HOME'), 'taylor_venv', 'bin', 'python'), 'python3', 'python'};
    py = '';
    for i = 1:numel(cands)
        c = cands{i}; if isempty(c), continue; end
        [st, ~] = system(sprintf('"%s" -c "import onnx, onnxruntime, vnnlib"', c));
        if st == 0, py = c; return; end
    end
end

function p = writeSat(xs)
    p = [tempname '.txt'];
    fid = fopen(p, 'w'); fprintf(fid, 'sat\n');
    for i = 1:numel(xs), fprintf(fid, '(X_%d %g)\n', i - 1, xs(i)); end
    fclose(fid);
end
