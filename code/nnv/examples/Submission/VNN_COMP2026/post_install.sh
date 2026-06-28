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

# ---- 2) Dedicated python venv (NNV_ORT_PYTHON) with matlab.engine + ONNX stack + torch ----
# The eval box has anaconda FIRST on PATH (smoke 269), so a bare `python3` is ambiguous and the R2026a
# MATLAB engine did NOT install into it. Build a self-contained venv from a SYSTEM python the engine
# supports, put EVERYTHING in it, and point the whole harness at it via NNV_ORT_PYTHON. Works anywhere.
echo "=== python diagnostics (eval box) ==="
echo "PATH python3 -> $(command -v python3) ($(python3 --version 2>&1))"
ls -d /home/ubuntu/anaconda3 2>/dev/null && echo "  (anaconda present on PATH)"
echo "/usr/bin/python3 -> $(/usr/bin/python3 --version 2>&1)"
for v in 3.13 3.12 3.11 3.10; do command -v python$v >/dev/null 2>&1 && echo "  python$v -> $(python$v --version 2>&1)"; done
sudo apt-get install -y python3-venv python3-pip || true
# Pick a SYSTEM python (NOT anaconda) that the R2026a MATLAB engine supports.
BASEPY=""
for cand in /usr/bin/python3.12 /usr/bin/python3.11 /usr/bin/python3.10 /usr/bin/python3 python3; do
    if command -v "$cand" >/dev/null 2>&1; then BASEPY="$cand"; break; fi
done
echo "Using BASEPY=$BASEPY ($("$BASEPY" --version 2>&1))"
VENV="$HOME/.nnv_venv"
"$BASEPY" -m venv "$VENV"
VPY="$VENV/bin/python"
"$VPY" -m pip install --upgrade pip wheel
# matlab.engine into the venv (capture output so a failure is visible, not silent)
( cd /usr/local/matlab/extern/engines/python && "$VPY" -m pip install . ) 2>&1 | tail -25 || echo "WARN: matlab.engine pip install . non-zero"
# the ONNX importer + witness-gate stack + torch into the venv
"$VPY" -m pip install -r "$REQ"
"$VPY" -m pip install torch
# FALLBACK: the PLATFORM's prep step runs onnx2nnv.py under a bare `python3` (may not see NNV_ORT_PYTHON),
# so also put the ONNX stack in the system python3 (best-effort; PEP-668 boxes need --break-system-packages).
pip3 install -r "$REQ" 2>/dev/null || python3 -m pip install -r "$REQ" --break-system-packages 2>/dev/null || true
# FATAL GATE from a NEUTRAL cwd (the old gate ran from the engine source dir -> false positive on `import matlab`).
cd /tmp && "$VPY" -c "import matlab.engine, numpy, scipy, onnx, onnxruntime, onnxsim, onnxoptimizer, vnnlib" \
    || { echo "ERROR: venv ($VPY, base=$BASEPY) missing matlab.engine or the onnx stack" >&2; "$VPY" --version >&2; exit 1; }
# Persist NNV_ORT_PYTHON for later shells (run_instance also exports it via vnncomp2026_env.sh).
echo "export NNV_ORT_PYTHON=\"$VENV/bin/python\"" >> ~/.bashrc
echo "export NNV_ORT_PYTHON=\"$VENV/bin/python\"" >> ~/.profile
export NNV_ORT_PYTHON="$VENV/bin/python"

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
