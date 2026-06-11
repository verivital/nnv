function [ok, available] = validate_witness_onnx(onnx_path, x_flat, Hs, out_tol)
%VALIDATE_WITNESS_ONNX  Definitive SAT-witness guard: replay through onnxruntime.
%   Replays the FLAT (ONNX-order) witness through the ORIGINAL ONNX model via
%   onnxruntime (tools/onnx2nnv_python/onnx_replay_check.py) and confirms the output
%   lands in the unsafe region G*y <= g for some Hs(h). This is the competition's
%   actual check, so unlike validate_witness (which replays through NNV's own
%   forward) it also catches a SYSTEMATIC import/encoding mismatch -- the likely
%   source of several of NNV's 19 incorrect/missing-CE (-150) penalties in 2025.
%   Use it as the final gate before emitting `sat`: downgrade to `unknown` ONLY when
%   `~ok && available` (onnxruntime ran and definitively found no violation); if
%   `~available` (onnxruntime unavailable / replay error), trust validate_witness.
%   See VNNCOMP2026_STRATEGY.md / VNNCOMP2026_PROGRESS.md (Pillar 2).
%
%   onnx_path : path to the ORIGINAL .onnx model.
%   x_flat    : witness input, FLAT in ONNX order (the value written to the CE file).
%   Hs        : array of HalfSpace (unsafe output region).
%   out_tol   : output-constraint slack (default 1e-4).
%
%   ok        : true only if onnxruntime confirms the witness violates some Hs(h).
%   available : true means the `ok` result is DEFINITIVE -- either a VIOLATED was found
%               (early return, ok=true), or EVERY Hs returned a clean OK (ok=false). So
%               `~ok && available` means onnxruntime definitively says the witness
%               violates NOTHING (the caller downgrades `sat`->`unknown`). available=false
%               means we could not check (Python/onnxruntime/model missing, or a replay
%               ERROR) -> the caller must TRUST validate_witness and NOT suppress the sat.

    if nargin < 4 || isempty(out_tol), out_tol = 1e-4; end
    ok = false; available = false;

    script = fullfile(nnvroot(), 'code', 'nnv', 'tools', 'onnx2nnv_python', 'onnx_replay_check.py');
    if ~isfile(script) || ~isfile(char(onnx_path))
        return;
    end

    td = tempname; mkdir(td);
    cleanup = onCleanup(@() rmdir(td, 's'));
    xf = fullfile(td, 'x_in.csv');
    writematrix(double(x_flat(:)), xf);

    % NB: file names must be case-DISTINCT -- Windows' filesystem is case-insensitive,
    % so "G.csv" and "g.csv" would be the same file and the second write clobbers the
    % first (g would overwrite G), corrupting the unsafe-region matrix.
    py = python_exe();
    all_clean = ~isempty(Hs);       % a "definitively non-violating" verdict needs >=1
    for h = 1:numel(Hs)             % Hs that onnxruntime cleanly checked
        Gf = fullfile(td, 'Gmat.csv'); gf = fullfile(td, 'gvec.csv');
        writematrix(double(Hs(h).G), Gf);
        writematrix(double(Hs(h).g(:)), gf);
        cmd = sprintf('%s "%s" "%s" "%s" "%s" "%s" --out-tol %g', ...
            py, script, char(onnx_path), xf, Gf, gf, out_tol);
        [st, out] = system(cmd); %#ok<ASGLU>
        if contains(out, 'VIOLATED')        % onnxruntime confirms a real counterexample
            ok = true; available = true; return;
        elseif ~contains(out, 'OK ')        % ERROR / python-missing -> cannot trust this
            all_clean = false;              % Hs (a clean non-violation prints "OK ...")
        end
    end
    available = all_clean;          % every Hs returned a clean OK -> confirmed non-CE
end

function py = python_exe()
    % Prefer an explicit interpreter if MATLAB knows one; else fall back to PATH.
    py = 'python';
    try
        pe = pyenv;
        if ~isempty(pe.Executable) && isfile(pe.Executable)
            py = ['"' char(pe.Executable) '"'];
        end
    catch
    end
end
