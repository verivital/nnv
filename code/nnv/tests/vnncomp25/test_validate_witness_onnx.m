function tests = test_validate_witness_onnx
%TEST_VALIDATE_WITNESS_ONNX  Tests for the onnxruntime SAT-witness guard.
%   validate_witness_onnx replays a FLAT (ONNX-order) witness through the ORIGINAL
%   .onnx via onnxruntime (the competition's own checker) and returns true ONLY when
%   the output lands in the unsafe region G*y <= g. It is the definitive Pillar-2
%   guard against the -150 "incorrect counterexample" penalty.
%
%   Two groups of tests:
%     * fallback contract (always runs): returns false -- never errors -- when the
%       model/Python is unavailable, so the runner safely degrades to `unknown`.
%     * replay correctness (runs only if python+onnxruntime are importable): exports
%       a tiny known-linear net to ONNX and checks violated->true / impossible->false.
    tests = functiontests(localfunctions);
end

% ---------- fallback contract (environment-independent) ----------

function test_missing_onnx_returns_false(tc)
    % A non-existent model path must yield false (safe), not an error.
    bogus = fullfile(tempdir, 'definitely_not_a_real_model_xyz.onnx');
    if isfile(bogus), delete(bogus); end
    [ok, available] = validate_witness_onnx(bogus, [0; 0], HalfSpace([1 0], 0));
    verifyFalse(tc, ok);
    verifyFalse(tc, available, 'missing model -> unavailable, so the runner trusts validate_witness (keeps sat)');
end

function test_empty_halfspace_array_returns_false(tc)
    % A VALID (existing) model path but an EMPTY unsafe region: the Hs loop never
    % iterates, so the wrapper returns false WITHOUT invoking Python. This exercises
    % the empty-region path itself (the missing-file guard is covered by the test
    % above). A stub file suffices -- no ONNX is parsed when there are zero HalfSpaces.
    f = [tempname '.onnx'];
    fid = fopen(f, 'w'); fwrite(fid, 'stub'); fclose(fid);
    c = onCleanup(@() delete(f));  %#ok<NASGU>
    ok = validate_witness_onnx(f, [1; 1], HalfSpace.empty);
    verifyFalse(tc, ok);
end

% ---------- replay correctness (needs python + onnxruntime) ----------

function test_violated_true_and_impossible_false(tc)
    assumeOnnxruntime(tc);
    [onnx, cleanup] = exportTinyLinearOnnx();  %#ok<ASGLU> keep cleanup alive
    % Net is y = 1*x1 + 2*x2 (single output). For x=[1;1], y=3.
    x = [1; 1];
    % Unsafe set {y : y <= 10} -> 3 <= 10 -> VIOLATED -> ok=true, available=true.
    [ok_v, av_v] = validate_witness_onnx(onnx, x, HalfSpace(1, 10));
    verifyTrue(tc, ok_v, 'witness y=3 should violate {y <= 10}');
    verifyTrue(tc, av_v, 'a clean VIOLATED verdict is available');
    % Unsafe set {y : y <= -10} -> impossible for y=3 -> OK -> ok=false, available=TRUE.
    % available=true + ok=false is the DOWNGRADE case: onnxruntime definitively says the
    % witness violates nothing, so the runner emits `unknown` instead of a wrong `sat`.
    [ok_i, av_i] = validate_witness_onnx(onnx, x, HalfSpace(1, -10));
    verifyFalse(tc, ok_i, 'witness y=3 must NOT be accepted for {y <= -10}');
    verifyTrue(tc, av_i, 'a clean OK (non-violation) verdict is available -> caller downgrades');
end

function test_multi_halfspace_any_violates(tc)
    assumeOnnxruntime(tc);
    [onnx, cleanup] = exportTinyLinearOnnx();  %#ok<ASGLU>
    x = [1; 1];   % y = 3
    Hs = [HalfSpace(1, -10), HalfSpace(1, 10)];  % 1st impossible, 2nd holds
    verifyTrue(tc, validate_witness_onnx(onnx, x, Hs), 'any-violates should accept');
end

% ---------- helpers ----------

function assumeOnnxruntime(tc)
    % Skip (not fail) the replay tests when onnxruntime is not importable, so CI
    % machines without it stay green. Resolve the interpreter the SAME way
    % validate_witness_onnx's python_exe() does (prefer MATLAB's configured pyenv,
    % else PATH 'python') so this gate matches the python the wrapper actually calls.
    [st, ~] = system([resolve_python() ' -c "import onnxruntime"']);
    assumeEqual(tc, st, 0, 'python/onnxruntime unavailable -- skipping replay tests');
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

function [onnxPath, cleanup] = exportTinyLinearOnnx()
    % A minimal known-linear net: featureInput(2) -> fc(1) with y = [1 2]*x.
    % Preset the weights ON THE LAYER before building the dlnetwork -- dlnetwork
    % preserves already-initialized learnables, so y is deterministic.
    fc = fullyConnectedLayer(1, 'Name', 'fc');
    fc.Weights = single([1 2]);   % y = 1*x1 + 2*x2
    fc.Bias    = single(0);
    layers = [featureInputLayer(2, 'Name', 'in'); fc];
    net = dlnetwork(layers);
    onnxPath = [tempname '.onnx'];
    exportONNXNetwork(net, onnxPath);
    cleanup = onCleanup(@() delete(onnxPath));
end
