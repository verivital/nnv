#!/bin/bash
# VNN-COMP 2026 post_install.sh for NNV. NOTE: the platform runs the post_install_tool *payload string*, NOT
# this repo file -- the submit payload MUST carry this file's content in d['post_install_tool'] (the empty
# field is exactly why smokes 266-277 did nothing: license, python, and MATLAB setup never ran). This file is
# the maintained SOURCE that the submitter uploads into that field.
#
# Runs AFTER install_tool.sh, on the SAME instance, BEFORE the per-benchmark runs (which execute as the
# default user `ubuntu`). With run_post_installation_script_as_root=True this runs as root (so the `sudo`
# calls below are harmless no-ops and license/apt/pip always have privilege). Order: DEBUG -> LICENSE ->
# system-python matlab.engine + ONNX stack -> licensed NNV install -> GPU hold.
set +e
# Mirror ALL output to a WORLD-WRITABLE /tmp log (any user can open it, unlike /home/ubuntu) so run_instance.sh
# (the one reliably-captured channel) can cat it back into the captured run log. Plain redirect, not the
# process-substitution that failed in smoke 273.
exec >> /tmp/nnv_post_install.log 2>&1
chmod 666 /tmp/nnv_post_install.log 2>/dev/null || true
echo "===================================================================================="
echo "=== POST_INSTALL START $(date -u +%FT%TZ) ==="
echo "whoami=$(whoami)  id=$(id)"
echo "HOME=$HOME  PWD=$(pwd)  SHELL=$SHELL"
echo "PATH=$PATH"

# ---- 0) EXHAUSTIVE ENVIRONMENT DEBUG (paths, MATLAB, python, sudo, license dir) ----
echo "--- [debug] sudo capability ---"
sudo -n true 2>&1; echo "sudo -n true rc=$?"
echo "--- [debug] MATLAB location ---"
echo "command -v matlab = $(command -v matlab 2>&1)"
echo "readlink matlab    = $(readlink -f "$(command -v matlab 2>/dev/null)" 2>&1)"
ls -ld /usr/local/matlab /usr/local/MATLAB /usr/local/MATLAB/* /opt/matlab* /opt/MATLAB* 2>&1 | sed 's/^/  /'
# Resolve the real MATLAB root robustly (used for license dir + engine source dir).
MLROOT="$(dirname "$(dirname "$(readlink -f "$(command -v matlab 2>/dev/null)" 2>/dev/null)")" 2>/dev/null)"
for cand in "$MLROOT" /usr/local/matlab /usr/local/MATLAB/R2026a /usr/local/MATLAB/R2026b /opt/matlab /opt/MATLAB/R2026a; do
    [ -n "$cand" ] && [ -d "$cand/extern/engines/python" ] && { MLROOT="$cand"; break; }
done
echo "MLROOT=$MLROOT"
ls -ld "$MLROOT" "$MLROOT/bin" "$MLROOT/licenses" "$MLROOT/extern/engines/python" 2>&1 | sed 's/^/  /'
echo "matlab -batch version:"; matlab -batch "disp(version); disp(matlabroot)" 2>&1 | head -5
echo "--- [debug] python interpreters ---"
echo "PATH python3 = $(command -v python3 2>&1) ($(python3 --version 2>&1))"
for p in /usr/bin/python3 /usr/bin/python3.10 /usr/bin/python3.11 /usr/bin/python3.12 /usr/bin/python3.13; do
    command -v "$p" >/dev/null 2>&1 && echo "  $p -> $("$p" --version 2>&1)"
done
echo "pip3 = $(command -v pip3 2>&1); uv = $(command -v uv 2>&1)"
echo "--- [debug] existing license dir BEFORE install ---"
ls -l "$MLROOT/licenses/" 2>&1 | sed 's/^/  /'

# Resolve repo root from THIS script's directory inside the clone (post_install_tool.sh is written into the
# scripts dir, so $0 resolves there); fall back to the known clone path.
SD="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd)"
[ -f "$SD/prepare_run.m" ] || SD="/home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2026"
NNV_ROOT="$(cd "$SD/../../.." 2>/dev/null && pwd)"
[ -f "$NNV_ROOT/install.m" ] || NNV_ROOT="/home/ubuntu/toolkit/code/nnv"
REQ="$NNV_ROOT/tools/onnx2nnv_python/requirements.txt"
echo "SD=$SD  NNV_ROOT=$NNV_ROOT  REQ=$REQ ($([ -f "$REQ" ] && echo present || echo MISSING))"

# ---- 1) LICENSE (first: the MATLAB install in step 3 needs a licensed MATLAB) ----
echo "=== [1] LICENSE ==="
mkdir -p ~/.matlab/R2026a_licenses
cd ~/.matlab/R2026a_licenses
# MAC-locked MATLAB R2026a license (HOSTID 02e1e896fadb == ENI eni-0b11771dfe21b94ee; exp 30-may-2027). URL QUOTED.
curl --retry 100 --retry-connrefused -fL -o license.lic "https://www.dropbox.com/scl/fi/w5jgddmf3qm5znjw67ajm/matlab-license-vnncomp2026-nnv.lic?rlkey=z3wnimbad4ykjyde95yhq7ik3&st=2oc64st3&dl=1"
echo "curl rc=$?"; sleep 3
ls -al license.lic 2>&1 | sed 's/^/  /'
echo "license head:"; head -3 license.lic 2>&1 | sed 's/^/  /'
if [ ! -s license.lic ] || ! grep -qE 'INCREMENT|MathWorks license|SERVER|DAEMON' license.lic; then
    echo "ERROR: MATLAB license download failed or invalid"
fi
sudo mkdir -p "$MLROOT/licenses"
sudo cp -f license.lic "$MLROOT/licenses/" && echo "license copied to $MLROOT/licenses/" || echo "ERROR: license cp failed"
sudo rm -f "$MLROOT/licenses/license_info.xml"
echo "license dir AFTER copy:"; sudo ls -l "$MLROOT/licenses/" 2>&1 | sed 's/^/  /'

# ---- 2) system-python matlab.engine + ONNX stack + torch (sudo --break-system-packages -> system site-packages) ----
echo "=== [2] PYTHON: matlab.engine + ONNX stack into a SYSTEM python ==="
sudo apt-get install -y python3-pip software-properties-common 2>&1 | tail -4
ENGSRC="$MLROOT/extern/engines/python"
echo "engine source dir: $ENGSRC ($([ -d "$ENGSRC" ] && echo present || echo MISSING))"
ls -l "$ENGSRC" 2>&1 | sed 's/^/  /'
REQ_OK() { ( cd /tmp && "$1" -c "import matlab.engine,numpy,scipy,onnx,onnxruntime,onnxsim,onnxoptimizer,vnnlib" ) >/dev/null 2>&1; }
WORKING_PY=""
for PY in /usr/bin/python3.13 /usr/bin/python3.12 /usr/bin/python3.11 /usr/bin/python3.10 /usr/bin/python3; do
    command -v "$PY" >/dev/null 2>&1 || continue
    echo "--- installing matlab.engine + deps into $PY ($("$PY" --version 2>&1)) ---"
    ( cd "$ENGSRC" && sudo "$PY" -m pip install . --break-system-packages ) 2>&1 | tail -20
    echo "  engine import test: $( cd /tmp && "$PY" -c 'import matlab.engine; print("ENGINE_OK")' 2>&1 | tail -1)"
    sudo "$PY" -m pip install -r "$REQ" --break-system-packages 2>&1 | tail -5
    sudo "$PY" -m pip install torch --break-system-packages 2>&1 | tail -3
    if REQ_OK "$PY"; then WORKING_PY="$PY"; echo "ENGINE_WORKS $PY"; break; fi
    echo "  full-stack import after install ($PY):"; ( cd /tmp && "$PY" -c "import matlab.engine,numpy,scipy,onnx,onnxruntime,onnxsim,onnxoptimizer,vnnlib" ) 2>&1 | tail -4 | sed 's/^/    /'
done
if [ -z "$WORKING_PY" ]; then
    echo "=== no pre-installed system python worked; installing python3.12 via deadsnakes ==="
    sudo add-apt-repository -y ppa:deadsnakes/ppa 2>&1 | tail -3
    sudo apt-get update 2>&1 | tail -3
    sudo apt-get install -y python3.12 python3.12-venv python3.12-dev python3.12-distutils 2>&1 | tail -4
    PY=/usr/bin/python3.12
    if command -v "$PY" >/dev/null 2>&1; then
        echo "deadsnakes $PY -> $("$PY" --version 2>&1)"
        sudo "$PY" -m ensurepip 2>&1 | tail -3
        ( cd "$ENGSRC" && sudo "$PY" -m pip install . --break-system-packages ) 2>&1 | tail -20
        sudo "$PY" -m pip install -r "$REQ" --break-system-packages 2>&1 | tail -5
        sudo "$PY" -m pip install torch --break-system-packages 2>&1 | tail -3
        REQ_OK "$PY" && WORKING_PY="$PY" && echo "ENGINE_WORKS $PY (deadsnakes)"
        [ -z "$WORKING_PY" ] && { echo "deadsnakes import diag:"; ( cd /tmp && "$PY" -c "import matlab.engine,numpy,onnx" ) 2>&1 | tail -4 | sed 's/^/    /'; }
    else
        echo "ERROR: deadsnakes python3.12 not present after install"
    fi
fi
# Best-effort: seed the platform-prep python3 (prepare_instance runs onnx2nnv.py under bare python3).
sudo python3 -m pip install -r "$REQ" --break-system-packages 2>/dev/null || python3 -m pip install -r "$REQ" 2>/dev/null || true
echo "WORKING_SYS_PY=${WORKING_PY:-NONE}"

# ---- 3) NNV install (tbxmanager: mpt/glpk/sedumi + savepath) -- NOW LICENSED ----
echo "=== [3] NNV install + prepare_run ==="
matlab -batch "cd('${NNV_ROOT}'); install" 2>&1 | tail -25 || echo "WARN: NNV tbxmanager install returned non-zero"
matlab -nodisplay -r "cd('${SD}'); prepare_run; quit" 2>&1 | tail -15 || echo "WARN: prepare_run returned non-zero"
# Ensure root-created NNV/toolbox files are world-readable for the run user `ubuntu`.
sudo chmod -R a+rX "$NNV_ROOT/tbxmanager" "$NNV_ROOT/code" 2>/dev/null || true

# ---- 4) GPU persistence + lock the driver ----
echo "=== [4] GPU ==="
nvidia-smi 2>&1 | head -5
sudo nvidia-smi -pm 1 2>/dev/null || true
NVPKG="$(dpkg -l 2>/dev/null | awk '/^ii +nvidia-driver-[0-9]/{print $2}' | head -1)"
sudo apt-mark hold linux-image-generic linux-headers-generic "${NVPKG:-nvidia-driver-570}" 2>/dev/null || true
sudo systemctl disable unattended-upgrades 2>/dev/null || true
echo "=== POST_INSTALL END $(date -u +%FT%TZ) WORKING_SYS_PY=${WORKING_PY:-NONE} MLROOT=$MLROOT ==="
