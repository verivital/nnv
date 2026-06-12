function tests = test_cctsdb_enumerate
%TEST_CCTSDB_ENUMERATE  Tests for the cctsdb_yolo complete-enumeration verifier.
%   cctsdb_enumerate.py exploits the cctsdb_yolo_2023 structure (only
%   X_12288/X_12289 free in [0,62]; the ONNX consumes them ONLY through
%   Cast(int64) truncation, so the output is piecewise-constant on unit cells)
%   to decide instances by enumerating the <= 3969 integer points through
%   onnxruntime: SOUND AND COMPLETE for this benchmark, UNKNOWN on any
%   structural deviation. Exit protocol: 10 = SAT (+witness CSV), 11 = UNSAT,
%   2 = UNKNOWN.
%
%   Three groups of tests:
%     * the python script exists next to the runner (always runs);
%     * known-instance replays (run only if python+onnxruntime are importable
%       AND the local vnncomp2026_benchmarks clone exists): patch-1 idx_00559
%       -> UNSAT (min Y_0 ~ 0.8128), patch-3 idx_00055 -> SAT at (0,23);
%     * the run_vnncomp_instance dispatcher no longer rejects cctsdb_yolo with
%       the old stub error "Working on supporting this one".
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    % The runner + script live in examples/Submission/VNN_COMP2025. The matrix CI
    % splits tests into shards and nothing guarantees that dir is on the path, so
    % add it here to be self-sufficient.
    here = fileparts(mfilename('fullpath'));
    sub = fullfile(here, '..', '..', 'examples', 'Submission', 'VNN_COMP2025');
    tc.TestData.subdir = sub;
    tc.TestData.addedPath = '';
    if isfolder(sub) && ~contains(path, sub)
        addpath(sub);
        tc.TestData.addedPath = sub;
    end
end

function teardownOnce(tc)
    if isfield(tc.TestData, 'addedPath') && ~isempty(tc.TestData.addedPath)
        rmpath(tc.TestData.addedPath);
    end
end

% ---------- (a) script presence (environment-independent) ----------

function test_script_exists(tc)
    verifyTrue(tc, isfile(fullfile(tc.TestData.subdir, 'cctsdb_enumerate.py')), ...
        'cctsdb_enumerate.py must live next to run_vnncomp_instance.m');
end

% ---------- (b) known-instance replays (needs python+onnxruntime+benchmarks) ----------

function test_known_unsat_instance(tc)
    % patch-1 idx_00559: the full 63x63 grid stays strictly ABOVE the unsafe
    % region Y_0 <= 0.5 (min Y_0 ~ 0.8128) -> a COMPLETE unsat proof, exit 11.
    [onnx, vnnlib] = locate_instance(tc, 'patch-1', 'spec_onnx_patch-1_idx_00559_0');
    witness = [tempname '.csv'];
    c = onCleanup(@() cleanup_file(witness));
    [st, out] = run_enumerate(tc, onnx, vnnlib, witness);
    verifyEqual(tc, st, 11, sprintf('expected UNSAT exit 11, got %d: %s', st, out));
    tok = regexp(out, 'UNSAT\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)', 'tokens', 'once');
    verifyNotEmpty(tc, tok, 'stdout must carry "UNSAT ymin ymax"');
    ymin = str2double(tok{1});
    verifyEqual(tc, ymin, 0.8128, 'AbsTol', 1e-3, ...
        'min Y_0 over the grid must reproduce the known 0.8128');
    verifyFalse(tc, isfile(witness), 'no witness CSV may be written for UNSAT');
end

function test_known_sat_instance(tc)
    % patch-3 idx_00055: cell (0,23) lands in the unsafe region (Y_0 = 0) ->
    % SAT with that point as witness, exit 10, full input vector CSV written.
    [onnx, vnnlib] = locate_instance(tc, 'patch-3', 'spec_onnx_patch-3_idx_00055_0');
    witness = [tempname '.csv'];
    c = onCleanup(@() cleanup_file(witness));
    [st, out] = run_enumerate(tc, onnx, vnnlib, witness);
    verifyEqual(tc, st, 10, sprintf('expected SAT exit 10, got %d: %s', st, out));
    tok = regexp(out, '\<SAT\s+(-?\d+)\s+(-?\d+)\s+([-+0-9.eE]+)', 'tokens', 'once');
    verifyNotEmpty(tc, tok, 'stdout must carry "SAT p q y"');
    p = str2double(tok{1}); q = str2double(tok{2}); y = str2double(tok{3});
    verifyEqual(tc, [p q], [0 23], 'witness cell must reproduce the known (0,23)');
    verifyLessThanOrEqual(tc, y, 0.5, 'witness output must be in the unsafe region');
    % witness CSV: FULL input vector, flat vnnlib order, free coords = (p,q)
    verifyTrue(tc, isfile(witness), 'SAT must write the witness CSV');
    x = readmatrix(witness);
    verifyEqual(tc, numel(x), 12296, 'witness must be the full 12296 input vector');
    verifyTrue(tc, all(isfinite(x)), 'witness values must all be finite');
    verifyEqual(tc, x(12288 + 1), p, 'X_12288 must equal the reported p'); % 1-based
    verifyEqual(tc, x(12289 + 1), q, 'X_12289 must equal the reported q');
end

% ---------- (c) dispatcher routing (mirrors classify_dispatch in the smoke test) ----------

function test_dispatcher_no_longer_stub_errors(tc)
    % run_vnncomp_instance must neither raise the old cctsdb stub error
    % ("Working on supporting this one") nor reject the category as
    % unsupported: it routes to the enumeration path, which degrades bogus
    % inputs to a sound 'unknown' WITHOUT raising.
    tc.assumeTrue(exist('run_vnncomp_instance', 'file') == 2, ...
        'run_vnncomp_instance not on path; add examples/Submission/VNN_COMP2025');
    onnx = [tempname '.onnx'];          % does not exist
    vnnlib = [tempname '.vnnlib'];      % does not exist
    outf = [tempname '.txt'];
    c = onCleanup(@() cleanup_file(outf));
    msg = '';
    try
        run_vnncomp_instance('cctsdb_yolo', onnx, vnnlib, outf);
    catch ME
        msg = ME.message;
    end
    verifyFalse(tc, contains(msg, 'Working on supporting this one'), ...
        'the cctsdb_yolo stub error must be gone');
    verifyFalse(tc, contains(msg, 'ONNX model not supported'), ...
        'cctsdb_yolo must not be rejected as an unsupported category');
    verifyEmpty(tc, msg, sprintf('dispatch must not raise at all; got: %s', msg));
    verifyTrue(tc, isfile(outf), 'runner must write the output file');
    verifyEqual(tc, strtrim(fileread(outf)), 'unknown', ...
        'bogus inputs must degrade to a sound ''unknown''');
end

% ---------- helpers ----------

function [st, out] = run_enumerate(tc, onnx, vnnlib, witness)
    % Invoke the script the same way verify_cctsdb_enumeration does. It accepts
    % .gz transparently, so the benchmark files need no decompression here.
    script = fullfile(tc.TestData.subdir, 'cctsdb_enumerate.py');
    cmd = sprintf('%s "%s" "%s" "%s" "%s" --timeout 300', ...
        resolve_python(), script, onnx, vnnlib, witness);
    [st, out] = system(cmd);
end

function [onnx, vnnlib] = locate_instance(tc, model, spec)
    % Skip (not fail) when python/onnxruntime or the local benchmark clone is
    % unavailable, so CI machines without them stay green.
    [pyok, ~] = system([resolve_python() ' -c "import onnx, onnxruntime, numpy"']);
    assumeEqual(tc, pyok, 0, 'python/onnx/onnxruntime unavailable -- skipping replay tests');
    here = fileparts(mfilename('fullpath'));
    root = getenv('NNV_VNNCOMP2026_BENCHMARKS');   % points at .../benchmarks
    if isempty(root)
        % default: sibling clone of the repo, <repo>/../vnncomp2026_benchmarks
        root = fullfile(here, '..', '..', '..', '..', '..', ...
            'vnncomp2026_benchmarks', 'benchmarks');
    end
    % accept the env var pointing at either the repo root or its benchmarks/ dir
    if ~isfolder(fullfile(root, 'cctsdb_yolo_2023')) && isfolder(fullfile(root, 'benchmarks', 'cctsdb_yolo_2023'))
        root = fullfile(root, 'benchmarks');
    end
    bench = fullfile(root, 'cctsdb_yolo_2023', '1.0');
    onnx = fullfile(bench, 'onnx', [model '.onnx.gz']);
    vnnlib = fullfile(bench, 'vnnlib', [spec '.vnnlib.gz']);
    assumeTrue(tc, isfile(onnx) && isfile(vnnlib), ...
        sprintf('benchmark clone not found under %s -- skipping replay tests', bench));
end

function py = resolve_python()
    % Same resolution idiom as validate_witness_onnx / verify_cctsdb_enumeration:
    % prefer MATLAB's configured pyenv executable, else PATH 'python'.
    py = 'python';
    try
        pe = pyenv;
        if ~isempty(pe.Executable) && isfile(pe.Executable)
            py = ['"' char(pe.Executable) '"'];
        end
    catch
    end
end

function cleanup_file(f)
    if exist(f, 'file'), delete(f); end
end
