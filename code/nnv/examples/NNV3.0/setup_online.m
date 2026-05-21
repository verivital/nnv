% setup_online.m -- one-time NNV setup for the nnv3.0-online image.
%
% Run this in the BROWSER MATLAB session at http://localhost:8888 after
% signing in to your MathWorks account. Configures the Python venv,
% installs NNV paths, and verifies AIVL availability. AIVL itself is
% already installed at build time by Dockerfile.online's `mpm install
% ... Deep_Learning_Toolbox_Verification_Library ...` -- no manual
% tarball staging or `toolbox_install.m` extraction needed.
%
% Usage (paste in the MATLAB Command Window):
%   run('/home/matlab/nnv/code/nnv/examples/NNV3.0/setup_online.m')
%
% Subsequent docker run invocations of nnv3.0-online with the same
% nnv3-matlab-prefs and nnv3-matlab-mw named volumes mounted reuse the
% saved MATLAB path. The script is idempotent: NNV install.m is
% safe to re-run.

NNV_ROOT = '/home/matlab/nnv/code/nnv';
VENV_PY  = '/home/matlab/nnv/.venv/bin/python';

fprintf('\n=== setup_online (1/2) Python venv ===\n');
try
    pyenv('Version', VENV_PY);
    pe = pyenv;
    fprintf('  Python: %s -> %s\n', pe.Version, pe.Executable);
catch ME
    fprintf(2, '  [pyenv warning] %s\n', ME.message);
end

fprintf('\n=== setup_online (2/2) NNV paths ===\n');
try
    cd(NNV_ROOT);
    install
    savepath
    cd(NNV_ROOT);
    check_nnv_setup
catch ME
    fprintf(2, '  [install warning] %s\n', ME.message);
end

fprintf('\n=== AIVL availability check ===\n');
vnr_path = which('verifyNetworkRobustness');
if ~isempty(vnr_path) && contains(vnr_path, 'aivnv')
    fprintf('  OK -- AIVL available at %s\n', vnr_path);
    fprintf('  ToolComparison will include MathWorks-side rows.\n');
else
    fprintf(2, '  AIVL NOT found on the MATLAB path.\n');
    fprintf(2, '  This means your MathWorks licence may not include the\n');
    fprintf(2, '  Deep Learning Toolbox Verification Library entitlement.\n');
    fprintf(2, '  ToolComparison will run NNV-only; MathWorks rows absent from table_main.\n');
end

fprintf('\n=== setup_online DONE ===\n');
fprintf('Next: run the smoke (~20-25 min) or full (~3-5 h) experiments.\n');
fprintf('  run(''/home/matlab/nnv/code/nnv/examples/NNV3.0/run_smoke.m'')\n');
fprintf('  run(''/home/matlab/nnv/code/nnv/examples/NNV3.0/run_full.m'')\n\n');
