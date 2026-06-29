#!/bin/bash
# VNN-COMP 2026 post_install.sh for NNV. Runs AFTER install_tool.sh, on the SAME instance, BEFORE the per-
# benchmark runs (which execute as the default user `ubuntu`). Does the heavy LICENSED setup, in order:
# LICENSE -> system-python matlab.engine + ONNX stack -> licensed NNV install -> GPU driver hold.
#
# TWO HARD LESSONS from smoke 266-276 (do NOT regress):
# 1) DO NOT redirect this script's output to /home/ubuntu/* . post_install runs as a user that CANNOT write
#    /home/ubuntu, so a `{ ...; } >> /home/ubuntu/log` redirect FAILS TO OPEN and the whole block then does
#    NOT execute -- that is why the earlier license/python/matlab setup silently never happened and no log
#    ever appeared. We now let output flow to the platform's captured stdout (no self-redirect).
# 2) DO NOT use a venv or any ~/.nnv_* artifact for the python: the runs are a DIFFERENT user (ubuntu) and
#    a venv symlinks to a per-user base python they cannot reach. Instead install matlab.engine + deps into a
#    SYSTEM python via `sudo --break-system-packages` -- system site-packages are root-owned, hence readable
#    and executable by EVERY user incl. ubuntu, and persist on the single shared instance. The runtime
#    (vnncomp2026_env.sh, run as ubuntu) then PROBES the system pythons and picks the one with matlab.engine.
set +e
echo "=== POST_INSTALL START $(date -u +%FT%TZ) whoami=$(whoami) HOME=$HOME pwd=$(pwd) ==="
id

# Resolve the cloned repo root from THIS script's location (robust to the clone path).
SD="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"          # .../code/nnv/examples/Submission/VNN_COMP2026
NNV_ROOT="$(cd "$SD/../../.." && pwd)"                          # .../code/nnv  (has install.m)
REQ="$NNV_ROOT/tools/onnx2nnv_python/requirements.txt"
echo "SD=$SD NNV_ROOT=$NNV_ROOT REQ=$REQ"

# ---- 1) LICENSE (must be first: the MATLAB install in step 3 needs a licensed MATLAB) ----
mkdir -p ~/.matlab/R2026a_licenses
cd ~/.matlab/R2026a_licenses
# MAC-locked MATLAB R2026a license (HOSTID 02e1e896fadb == ENI eni-0b11771dfe21b94ee; exp 30-may-2027).
# URL QUOTED so the &-containing query string is not split / backgrounded.
curl --retry 100 --retry-connrefused -L -o license.lic "https://www.dropbox.com/scl/fi/w5jgddmf3qm5znjw67ajm/matlab-license-vnncomp2026-nnv.lic?rlkey=z3wnimbad4ykjyde95yhq7ik3&st=2oc64st3&dl=1"
sleep 5
ls -al
if [ ! -s license.lic ] || ! grep -qE 'INCREMENT|MathWorks license' license.lic; then
    echo "ERROR: MATLAB license download failed or invalid (missing/empty/not a passcode file)"
fi
sudo cp -f license.lic /usr/local/matlab/licenses/ && echo "license installed to /usr/local/matlab/licenses/" || echo "ERROR: failed to install license"
[ -s /usr/local/matlab/licenses/license.lic ] && echo "license present after copy" || echo "ERROR: license not present after copy"
sudo rm -f /usr/local/matlab/licenses/license_info.xml

# ---- 2) system-python matlab.engine + ONNX stack + torch (sudo --break-system-packages -> system site-packages) ----
echo "=== system-python engine install ==="
echo "PATH python3 -> $(command -v python3) ($(python3 --version 2>&1)); /usr/bin/python3 -> $(/usr/bin/python3 --version 2>&1)"
sudo apt-get install -y python3-pip software-properties-common 2>&1 | tail -3
REQ_OK() { ( cd /tmp && "$1" -c "import matlab.engine,numpy,scipy,onnx,onnxruntime,onnxsim,onnxoptimizer,vnnlib" ) >/dev/null 2>&1; }
WORKING_PY=""
for PY in /usr/bin/python3.13 /usr/bin/python3.12 /usr/bin/python3.11 /usr/bin/python3.10 /usr/bin/python3; do
    command -v "$PY" >/dev/null 2>&1 || continue
    echo "--- installing matlab.engine + deps into $PY ($("$PY" --version 2>&1)) ---"
    ( cd /usr/local/matlab/extern/engines/python && sudo "$PY" -m pip install . --break-system-packages ) 2>&1 | tail -15
    sudo "$PY" -m pip install -r "$REQ" --break-system-packages 2>&1 | tail -4
    sudo "$PY" -m pip install torch --break-system-packages 2>&1 | tail -2
    if REQ_OK "$PY"; then WORKING_PY="$PY"; echo "ENGINE_WORKS $PY"; break; fi
done
# If no pre-installed system python supports the R2026a engine, add python3.12 system-wide via deadsnakes + retry.
if [ -z "$WORKING_PY" ]; then
    echo "=== no pre-installed system python worked; installing python3.12 via deadsnakes ==="
    sudo add-apt-repository -y ppa:deadsnakes/ppa 2>&1 | tail -2
    sudo apt-get update 2>&1 | tail -2
    sudo apt-get install -y python3.12 python3.12-venv python3.12-dev python3.12-distutils 2>&1 | tail -3
    PY=/usr/bin/python3.12
    if command -v "$PY" >/dev/null 2>&1; then
        sudo "$PY" -m ensurepip 2>&1 | tail -2
        ( cd /usr/local/matlab/extern/engines/python && sudo "$PY" -m pip install . --break-system-packages ) 2>&1 | tail -15
        sudo "$PY" -m pip install -r "$REQ" --break-system-packages 2>&1 | tail -4
        sudo "$PY" -m pip install torch --break-system-packages 2>&1 | tail -2
        REQ_OK "$PY" && WORKING_PY="$PY" && echo "ENGINE_WORKS $PY (deadsnakes)"
    fi
fi
# Best-effort: also seed the platform-prep python3 (prepare_instance runs onnx2nnv.py under bare python3).
sudo python3 -m pip install -r "$REQ" --break-system-packages 2>/dev/null || python3 -m pip install -r "$REQ" 2>/dev/null || true
echo "WORKING_SYS_PY=${WORKING_PY:-NONE}"

# ---- 3) NNV install (tbxmanager: mpt/glpk/sedumi + savepath) -- NOW LICENSED ----
matlab -batch "cd('${NNV_ROOT}'); install" || echo "WARN: NNV tbxmanager install returned non-zero"
matlab -nodisplay -r "cd('${SD}'); prepare_run; quit" || echo "WARN: prepare_run returned non-zero"

# ---- 4) GPU persistence + lock the driver (pair with the form's restart-after-post-install) ----
sudo nvidia-smi -pm 1 2>/dev/null || true
NVPKG="$(dpkg -l 2>/dev/null | awk '/^ii +nvidia-driver-[0-9]/{print $2}' | head -1)"
sudo apt-mark hold linux-image-generic linux-headers-generic "${NVPKG:-nvidia-driver-570}" 2>/dev/null || true
sudo systemctl disable unattended-upgrades 2>/dev/null || true
echo "=== POST_INSTALL END $(date -u +%FT%TZ) WORKING_SYS_PY=${WORKING_PY:-NONE} ==="
