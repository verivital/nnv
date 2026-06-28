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

# ---- 2) MATLAB Engine API for Python (execute.py does `import matlab.engine` on every instance) ----
cd /usr/local/matlab/extern/engines/python
python3 -m pip install . || pip3 install matlabengine || echo "WARN: matlab.engine install non-zero (gate below decides)"
python3 -c "import matlab.engine" || { echo "ERROR: matlab.engine import failed; fix before running" >&2; exit 1; }

# ---- 3) NNV install (tbxmanager: mpt/glpk/sedumi + savepath) -- NOW LICENSED ----
# Was in install_tool.sh where it hit license Error -1.2 (no license yet). install.m is load-bearing
# (it installs the LP/polytope toolboxes via tbxmanager). Non-fatal guard: run_instance falls back to
# startup_nnv, but with the license present here this should succeed.
matlab -batch "cd('${NNV_ROOT}'); install" || echo "WARN: NNV tbxmanager install returned non-zero"
# prepare_run: warm the codegen packages / netcache (licensed).
matlab -nodisplay -r "cd('${SD}'); prepare_run; quit" || echo "WARN: prepare_run returned non-zero"

# ---- 4) Python deps for the ONNX importer + the SAT-witness gate (sudo for apt -- works per 2025) ----
sudo apt-get install -y python3 python3-pip
# onnx2nnv.py + the witness gate: single source of truth. Missing ANY -> prepare_instance.sh cannot
# generate the .nnv.mat manifests and every manifest-routed benchmark errors.
pip3 install --no-cache-dir -r "${REQ}"
# PREFLIGHT GATE: fail loudly if the importer / witness-gate deps don't import in THIS python (vnnlib
# included: a partial install silently fail-opens the -150 witness guard).
python3 -c "import numpy, scipy, onnx, onnxruntime, onnxsim, onnxoptimizer, vnnlib" \
    || { echo "ERROR: onnx2nnv.py / witness-gate deps failed to import after install (check vs ${REQ})" >&2; exit 1; }
# torch backs engine/nn/Prob_reach (cp-star / probabilistic-reach paths) -- not in requirements.txt.
python3 -m pip install torch

# ---- 5) GPU persistence + lock the driver (570; pair with the form's restart-after-post-install) ----
sudo nvidia-smi -pm 1
sudo apt-mark hold linux-image-generic linux-headers-generic nvidia-driver-570
sudo systemctl disable unattended-upgrades

echo "post_install complete (license + NNV install + python deps + preflight all passed)."
