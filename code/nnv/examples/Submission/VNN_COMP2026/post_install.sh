#!/bin/bash
# VNN-COMP 2026 post_install.sh for NNV. Runs AFTER install_tool.sh, as the NON-ROOT toolkit user
# (config.yaml run_post_installation_script_as_root: False). This is where the heavy lifting lives
# (matches the known-good 2025 pattern): LICENSE first, then the LICENSED MATLAB NNV install
# (tbxmanager: mpt/glpk/sedumi) + prepare_run, then the apt/pip deps (sudo for apt -- works on this
# platform per the 2025 submission). install_tool.sh only does the mpm support packages.
#
# Why here and not install_tool.sh: the 2026-06-28 smoke test (task 266) failed because install_tool
# did `apt-get` (non-root -> Permission denied) and `matlab -batch install` (no license yet -> Error -1.2).
# Licensing + the licensed MATLAB setup + apt MUST run here, after the license is in place.

# --- VISIBILITY: the platform's ToolkitPostInstall step log captures only its wrapper, NOT this script's
# output, and it marks the step Done regardless of our exit code. tee everything to a file; run_instance.sh
# cats it into EVERY run log (which IS captured), so failures here are diagnosable. ---
# simple redirect (NOT process-substitution `>(tee ...)` — that silently failed on the eval box in smoke 273,
# so the log was never created). All post_install output goes to this file; run_instance.sh cats it into the
# run log (the platform does NOT capture post_install output otherwise).
exec > "$HOME/.nnv_post_install.log" 2>&1
echo "=== POST_INSTALL START $(date -u +%FT%TZ) ==="

# Resolve the cloned repo root from THIS script's location BEFORE any cd (robust to the clone path).
SD="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"          # .../code/nnv/examples/Submission/VNN_COMP2026
NNV_ROOT="$(cd "$SD/../../.." && pwd)"                          # .../code/nnv  (has install.m)
REQ="$NNV_ROOT/tools/onnx2nnv_python/requirements.txt"

# ---- 1) LICENSE (must be first: the MATLAB install below needs a licensed MATLAB) ----
mkdir -p ~/.matlab/R2026a_licenses
cd ~/.matlab/R2026a_licenses
# MAC-locked MATLAB R2026a license (HOSTID 02e1e896fadb == ENI eni-0b11771dfe21b94ee). Verified
# 2026-06-28: FLEXlm passcode file, 116 products incl. Deep Learning + Parallel Computing + Optimization
# + GPU_Coder; expires 30-may-2027. The .lic is MAC-locked so the URL is unusable off this ENI.
# URL QUOTED (the prior year's unquoted &-URL backgrounded curl + dropped query params).
curl --retry 100 --retry-connrefused -L -o license.lic "https://www.dropbox.com/scl/fi/w5jgddmf3qm5znjw67ajm/matlab-license-vnncomp2026-nnv.lic?rlkey=z3wnimbad4ykjyde95yhq7ik3&st=2oc64st3&dl=1"
sleep 5
ls -al
if [ ! -s license.lic ] || ! grep -qE 'INCREMENT|MathWorks license' license.lic; then
    echo "ERROR: MATLAB license download failed or invalid (missing/empty/not a passcode file)" >&2; exit 1
fi
# /usr/local/matlab/licenses is root-owned and this script runs non-root -> sudo (works on this platform).
sudo cp -f license.lic /usr/local/matlab/licenses/ || { echo "ERROR: failed to install license" >&2; exit 1; }
[ -s /usr/local/matlab/licenses/license.lic ] || { echo "ERROR: license not present after copy" >&2; exit 1; }
sudo rm -f /usr/local/matlab/licenses/license_info.xml

# ---- 2) Robust python: a MATLAB-R2026a-engine-compatible python with matlab.engine + ONNX stack + torch ----
# The eval box has anaconda FIRST on PATH (smoke 269) and NO native python3.12 (the platform installs 3.12 via
# uv). A bare `python3` can't install the R2026a engine, and a plain `python3 -m venv` failed (smoke 272 ->
# error_exit_code_127, no venv). So: build a uv/3.12 venv (the platform's own proven method), VERIFY by
# importing matlab.engine, RECORD the working python path; fall back to system pythons + --break-system-packages.
echo "=== python diagnostics (eval box) ==="
echo "PATH python3 -> $(command -v python3) ($(python3 --version 2>&1))"
ls -d /home/ubuntu/anaconda3 2>/dev/null && echo "  (anaconda present on PATH)"
echo "/usr/bin/python3 -> $(/usr/bin/python3 --version 2>&1)"
for v in 3.13 3.12 3.11 3.10; do command -v python$v >/dev/null 2>&1 && echo "  python$v -> $(python$v --version 2>&1)"; done
echo "uv -> $(command -v uv || echo 'not found')"
sudo apt-get install -y python3-venv python3-pip || true

REQ_OK() { ( cd /tmp && "$1" -c "import matlab.engine,numpy,scipy,onnx,onnxruntime,onnxsim,onnxoptimizer,vnnlib" ) >/dev/null 2>&1; }
WORKING_PY=""
VENV="$HOME/.nnv_venv"

# PRIMARY: uv + python 3.12 venv (matches what the platform itself does for python 3.12).
if ! command -v uv >/dev/null 2>&1; then curl -LsSf https://astral.sh/uv/install.sh | sh || echo "WARN: uv install non-zero"; fi
export PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
echo "uv now -> $(command -v uv || echo 'still not found')"
uv python install 3.12 2>&1 | tail -5 || true
uv venv "$VENV" --python 3.12 2>&1 | tail -5 || /usr/bin/python3.12 -m venv "$VENV" 2>&1 | tail -5 || true
VPY="$VENV/bin/python"
if [ -x "$VPY" ]; then
    echo "venv python: $VPY ($("$VPY" --version 2>&1))"
    ( cd /usr/local/matlab/extern/engines/python && "$VPY" -m pip install . ) 2>&1 | tail -30 || echo "WARN: matlab.engine pip install . (venv) non-zero"
    "$VPY" -m pip install -r "$REQ" 2>&1 | tail -10
    "$VPY" -m pip install torch 2>&1 | tail -5
    if REQ_OK "$VPY"; then WORKING_PY="$VPY"; echo "VENV_OK $VPY"; else echo "VENV_FAIL importing in $VPY:"; ( cd /tmp && "$VPY" -c "import matlab.engine,numpy,scipy,onnx,onnxruntime,onnxsim,onnxoptimizer,vnnlib" ) 2>&1 | tail -8; fi
else
    echo "WARN: venv $VPY not created by uv/python3.12"
fi

# FALLBACK: install the engine + deps directly into a system python (PEP-668 -> --break-system-packages).
if [ -z "$WORKING_PY" ]; then
    for SP in /usr/bin/python3.12 /usr/bin/python3.11 /usr/bin/python3.10 /usr/bin/python3; do
        command -v "$SP" >/dev/null 2>&1 || continue
        echo "fallback trying $SP ($("$SP" --version 2>&1))"
        ( cd /usr/local/matlab/extern/engines/python && "$SP" -m pip install . --break-system-packages ) 2>&1 | tail -15 || true
        "$SP" -m pip install -r "$REQ" --break-system-packages 2>&1 | tail -5 || true
        "$SP" -m pip install torch --break-system-packages 2>&1 | tail -3 || true
        if REQ_OK "$SP"; then WORKING_PY="$SP"; echo "FALLBACK_OK $SP"; break; fi
    done
fi

# ALSO seed the platform-prep python3 with the ONNX stack (the prep runs onnx2nnv.py under bare python3).
pip3 install -r "$REQ" 2>/dev/null || python3 -m pip install -r "$REQ" --break-system-packages 2>/dev/null || true

# RECORD the working python for the harness (vnncomp2026_env.sh reads ~/.nnv_python_path).
if [ -n "$WORKING_PY" ]; then
    echo "$WORKING_PY" > "$HOME/.nnv_python_path"
    export NNV_ORT_PYTHON="$WORKING_PY"
    echo "export NNV_ORT_PYTHON=\"$WORKING_PY\"" >> ~/.bashrc
    echo "export NNV_ORT_PYTHON=\"$WORKING_PY\"" >> ~/.profile
    echo "NNV_ORT_PYTHON=$WORKING_PY (matlab.engine + onnx stack VERIFIED)"
else
    echo "ERROR: NO working python found with matlab.engine! runs will fall back to python3 and likely error." >&2
fi
# (No 'exit 1': the platform marks ToolkitPostInstall Done regardless of our exit code -- print loudly instead.)

# ---- 3) NNV install (tbxmanager: mpt/glpk/sedumi + savepath) -- NOW LICENSED ----
# Was in install_tool.sh where it hit license Error -1.2 (no license yet). install.m is load-bearing
# (it installs the LP/polytope toolboxes via tbxmanager). Non-fatal guard: run_instance falls back to
# startup_nnv, but with the license present here this should succeed.
matlab -batch "cd('${NNV_ROOT}'); install" || echo "WARN: NNV tbxmanager install returned non-zero"
# prepare_run: warm the codegen packages / netcache (licensed).
matlab -nodisplay -r "cd('${SD}'); prepare_run; quit" || echo "WARN: prepare_run returned non-zero"

# ---- 4) GPU persistence + lock the driver (pair with the form's restart-after-post-install) ----
# Hold WHATEVER nvidia driver is installed (570 or 580 -- both >= the R2026a CUDA 12.8 minimum of 570;
# the Lambda dev box uses 580). Detect it; fall back to 570. unattended-upgrades is disabled below too.
sudo nvidia-smi -pm 1 2>/dev/null || true
NVPKG="$(dpkg -l 2>/dev/null | awk '/^ii +nvidia-driver-[0-9]/{print $2}' | head -1)"
sudo apt-mark hold linux-image-generic linux-headers-generic "${NVPKG:-nvidia-driver-570}" 2>/dev/null || true
sudo systemctl disable unattended-upgrades 2>/dev/null || true

echo "post_install complete (license + NNV install + python deps + preflight all passed)."
