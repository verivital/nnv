% setup_online.m -- one-time NNV + AIVL setup for the nnv3.0-online image.
%
% Run this in the BROWSER MATLAB session at http://localhost:8888 after
% signing in to your MathWorks account. Equivalent to what install.m and
% toolbox_install.m do at build time on the network-licence Dockerfile,
% but deferred to runtime because matlab -batch cannot use a browser
% sign-in.
%
% Usage (paste in the MATLAB Command Window):
%   run('/home/matlab/nnv/code/nnv/examples/NNV3.0/setup_online.m')
%
% Subsequent docker run invocations of nnv3.0-online with the same
% nnv3-matlab-prefs and nnv3-matlab-mw named volumes mounted reuse the
% saved MATLAB path. AIVL extraction is idempotent: it skips if the
% support package is already present.

NNV_ROOT = '/home/matlab/nnv/code/nnv';
VENV_PY  = '/home/matlab/nnv/.venv/bin/python';
AIVL_INSTALLED_MARKER = '/home/matlab/Documents/MATLAB/SupportPackages/R2025b/toolbox/nnet/supportpackages/aivnv/verifyNetworkRobustness.m';

fprintf('\n=== setup_online (1/3) Python venv ===\n');
try
    pyenv('Version', VENV_PY);
    pe = pyenv;
    fprintf('  Python: %s -> %s\n', pe.Version, pe.Executable);
catch ME
    fprintf(2, '  [pyenv warning] %s\n', ME.message);
end

fprintf('\n=== setup_online (2/3) NNV paths ===\n');
try
    cd(NNV_ROOT);
    install
    savepath
    cd(NNV_ROOT);
    check_nnv_setup
catch ME
    fprintf(2, '  [install warning] %s\n', ME.message);
end

fprintf('\n=== setup_online (3/3) AIVL Support Package ===\n');
if exist(AIVL_INSTALLED_MARKER, 'file')
    fprintf('  AIVL already installed at %s -- skipping.\n', AIVL_INSTALLED_MARKER);
else
    try
        cd('/home/matlab/nnv/code/nnv/examples/NNV3.0/ToolComparison/utils');
        run('toolbox_install.m');
    catch ME
        fprintf(2, '  [toolbox_install warning] %s\n', ME.message);
    end
end

fprintf('\n=== AIVL availability check ===\n');
if exist(AIVL_INSTALLED_MARKER, 'file')
    fprintf('  OK -- ToolComparison will include MathWorks-side rows.\n');
else
    fprintf(2, '  AIVL NOT present at %s\n', AIVL_INSTALLED_MARKER);
    fprintf(2, '  ToolComparison will run NNV-only; MathWorks rows absent from table_main.\n');
    fprintf(2, '  Check that ToolComparison/utils/atva26-aivl.tar.gz exists in the image.\n');
end

fprintf('\n=== setup_online DONE ===\n');
fprintf('Next: run the smoke (~20-25 min) or full (~3-5 h) experiments.\n');
fprintf('  run(''/home/matlab/nnv/code/nnv/examples/NNV3.0/run_smoke.m'')\n');
fprintf('  run(''/home/matlab/nnv/code/nnv/examples/NNV3.0/run_full.m'')\n\n');
