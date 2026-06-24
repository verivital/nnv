function [verdict, detail] = authoritative_witness_gate(onnx, vnnlib, resultfile)
%AUTHORITATIVE_WITNESS_GATE  Re-validate an emitted `sat` result file against the AUTHORITATIVE
%   vnnlib parser + a real-onnx onnxruntime forward -- the SAME check execute.py runs on the
%   official path (validate_witness_authoritative.py). This is the run_all_benchmarks (dev-sweep)
%   twin of that gate: the sweep bypasses execute.py, so without this a sweep scorecard's `sat`
%   verdicts are only consensus-level-checked, never per-witness-checked. Wiring it into the sweep
%   makes a scorecard per-witness trustworthy (the 2026-06-24 sweep's 0-wrong was consensus-level).
%
%   verdict (char):
%     'valid'    -- the witness is a real counterexample on the real ONNX -> keep `sat`.
%     'spurious' -- some assertion fails on the real ONNX -> the caller MUST downgrade sat->unknown
%                   (a would-be -150 turned into 0 points). The ONLY non-fail-open outcome.
%     'cant'     -- could not check (missing deps / parse fail / shape mismatch / non-sat file /
%                   timeout) -> FAIL OPEN: the caller keeps today's verdict. The gate is then
%                   monotonically no-worse than not running it.
%
%   SOUND DIRECTION: the gate can only ever turn sat->unknown (lose +10), never manufacture a
%   verdict. The single source of truth for the actual check is validate_witness_authoritative.py
%   (exit 0=valid / 1=spurious / 2=cant); this wrapper only resolves a deps-complete python, bounds
%   the wall time, and maps the exit code. Disable via NNV_SWEEP_WITNESS_GATE=0 (caller's choice).
    verdict = 'cant'; detail = '';
    here = fileparts(mfilename('fullpath'));
    gate = fullfile(here, 'validate_witness_authoritative.py');
    if ~isfile(gate),        detail = 'validate_witness_authoritative.py missing'; return; end
    if ~isfile(char(resultfile)), detail = 'result file missing'; return; end

    py = i_ort_vnnlib_python();   % needs onnx+onnxruntime+vnnlib; '' if none on this box
    if isempty(py), detail = 'no python imports onnx+onnxruntime+vnnlib'; return; end

    % Bound the wall time so a hung onnxruntime can never stall the sweep loop. `timeout` exists on
    % the Linux eval/dev boxes; skip the prefix on Windows (rc 124 from timeout -> 'cant' = fail open).
    if ispc
        cmd = sprintf('%s "%s" "%s" "%s" "%s"', py, gate, char(onnx), char(vnnlib), char(resultfile));
    else
        cmd = sprintf('timeout 90 %s "%s" "%s" "%s" "%s"', py, gate, char(onnx), char(vnnlib), char(resultfile));
    end
    [rc, out] = system(cmd);
    detail = strtrim(out);
    switch rc
        case 0, verdict = 'valid';
        case 1, verdict = 'spurious';
        otherwise, verdict = 'cant';   % rc==2 (cant-check) or 124 (timeout) or anything else -> fail open
    end
end

function py = i_ort_vnnlib_python()
%I_ORT_VNNLIB_PYTHON  First interpreter that imports onnx+onnxruntime+vnnlib (quoted), else ''.
%   Mirrors run_vnncomp_instance/python_exe's probe order (NNV_ORT_PYTHON -> taylor_venv -> pyenv ->
%   system) but REQUIRES vnnlib too (the authoritative gate's parser). Kept SEPARATE from python_exe
%   on purpose: python_exe is shared with cctsdb_enumerate.py, which needs only onnx+ort -- requiring
%   vnnlib there would let a no-vnnlib box fall through to a no-ort interpreter (the cctsdb 39->0
%   regression). Cached for the MATLAB-process lifetime.
    persistent cached
    if ~isempty(cached), py = cached{1}; return; end
    cands = {};
    e = strtrim(getenv('NNV_ORT_PYTHON'));
    if ~isempty(e), cands{end+1} = e; end
    home = getenv('HOME');
    if ~isempty(home), cands{end+1} = fullfile(home, 'taylor_venv', 'bin', 'python'); end
    try
        pe = pyenv;
        if ~isempty(pe.Executable) && isfile(pe.Executable), cands{end+1} = char(pe.Executable); end
    catch
    end
    cands = [cands, {'python3', 'python'}];
    py = '';
    for i = 1:numel(cands)
        c = cands{i};
        if isempty(c), continue; end
        [st, ~] = system(sprintf('"%s" -c "import onnx, onnxruntime, vnnlib"', c));
        if st == 0, py = ['"' c '"']; break; end
    end
    cached = {py};
end
