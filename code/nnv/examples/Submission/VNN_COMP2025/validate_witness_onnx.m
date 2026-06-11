function ok = validate_witness_onnx(onnx_path, x_flat, Hs, out_tol)
%VALIDATE_WITNESS_ONNX  Definitive SAT-witness guard: replay through onnxruntime.
%   Replays the FLAT (ONNX-order) witness through the ORIGINAL ONNX model via
%   onnxruntime (tools/onnx2nnv_python/onnx_replay_check.py) and confirms the output
%   lands in the unsafe region G*y <= g for some Hs(h). This is the competition's
%   actual check, so unlike validate_witness (which replays through NNV's own
%   forward) it also catches a SYSTEMATIC import/encoding mismatch -- the likely
%   source of several of NNV's 19 incorrect/missing-CE (-150) penalties in 2025.
%   Use it as the final gate before emitting `sat`; emit `unknown` if it returns false.
%   See VNNCOMP2026_STRATEGY.md / VNNCOMP2026_PROGRESS.md (Pillar 2).
%
%   onnx_path : path to the ORIGINAL .onnx model.
%   x_flat    : witness input, FLAT in ONNX order (the value written to the CE file).
%   Hs        : array of HalfSpace (unsafe output region).
%   out_tol   : output-constraint slack (default 1e-4).
%
%   ok : true only if onnxruntime confirms the witness violates some Hs(h).
%   Returns false (safe) if Python/onnxruntime is unavailable or errors.

    if nargin < 4 || isempty(out_tol), out_tol = 1e-4; end
    ok = false;

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
    for h = 1:numel(Hs)
        Gf = fullfile(td, 'Gmat.csv'); gf = fullfile(td, 'gvec.csv');
        writematrix(double(Hs(h).G), Gf);
        writematrix(double(Hs(h).g(:)), gf);
        cmd = sprintf('%s "%s" "%s" "%s" "%s" "%s" --out-tol %g', ...
            py, script, char(onnx_path), xf, Gf, gf, out_tol);
        [st, out] = system(cmd);
        if st == 0 && contains(out, 'VIOLATED')
            ok = true; return;
        end
    end
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
