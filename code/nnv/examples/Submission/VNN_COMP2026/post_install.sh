#!/bin/bash
# VNN-COMP 2026 post_install.sh for NNV. NOTE: the platform runs the post_install_tool *payload string*, NOT
# this repo file -- the submit payload MUST carry this file's content in d['post_install_tool'] (the empty
# field is exactly why smokes 266-277 did nothing). This file is the maintained SOURCE the submitter uploads.
#
# Runs AFTER install_tool.sh, on the SAME instance, BEFORE the per-benchmark runs (which execute as user
# `ubuntu`). With run_post_installation_script_as_root=True this runs as ROOT (the `sudo` calls are harmless
# no-ops; license cp + system-wide pip always have privilege). Order: DEBUG -> LICENSE -> system-python
# matlab.engine + ONNX stack -> licensed NNV install -> GPU hold.
#
# VISIBILITY: do NOT redirect the body (no `exec >>file`). The platform DOES capture un-redirected body
# stdout/stderr in the ToolkitPostInstall step log (smoke 278 proved it: `+ set +e`/`+ exec` showed up), so
# letting output flow makes the whole post-install STREAM LIVE into that step log. The smoke-278 `exec >>/tmp`
# redirect turned post-install into a 35-min black box -- removed.
# HANG SAFETY: every `matlab` call and the license curl are wrapped in `timeout` so a stuck step fails fast
# (bounded) instead of hanging the whole task. SPEED: install the engine FIRST and skip the heavy onnx+torch
# install on any python whose matlab.engine does not import (avoids wasting a ~4-min torch install on a
# doomed interpreter); dedupe interpreters by realpath.
set +e
echo "===================================================================================="
echo "=== POST_INSTALL START $(date -u +%FT%TZ) ==="
echo "whoami=$(whoami)  id=$(id)"
echo "HOME=$HOME  PWD=$(pwd)  PATH=$PATH"

# ---- 0) ENVIRONMENT DEBUG (paths, MATLAB, python, sudo, license dir) ----
echo "--- [debug] sudo capability ---"; sudo -n true 2>&1; echo "sudo -n true rc=$?"
echo "--- [debug] MATLAB location ---"
echo "command -v matlab = $(command -v matlab 2>&1)"
echo "readlink matlab    = $(readlink -f "$(command -v matlab 2>/dev/null)" 2>&1)"
ls -ld /usr/local/matlab /usr/local/MATLAB /usr/local/MATLAB/* /opt/matlab* /opt/MATLAB* 2>&1 | sed 's/^/  /'
MLROOT="$(dirname "$(dirname "$(readlink -f "$(command -v matlab 2>/dev/null)" 2>/dev/null)")" 2>/dev/null)"
for cand in "$MLROOT" /usr/local/matlab /usr/local/MATLAB/R2026a /usr/local/MATLAB/R2026b /opt/matlab /opt/MATLAB/R2026a; do
    [ -n "$cand" ] && [ -d "$cand/extern/engines/python" ] && { MLROOT="$cand"; break; }
done
echo "MLROOT=$MLROOT"
ls -ld "$MLROOT" "$MLROOT/bin" "$MLROOT/licenses" "$MLROOT/extern/engines/python" 2>&1 | sed 's/^/  /'
echo "--- [debug] python interpreters ---"
echo "PATH python3 = $(command -v python3 2>&1) ($(python3 --version 2>&1))"
for p in /usr/bin/python3 /usr/bin/python3.10 /usr/bin/python3.11 /usr/bin/python3.12 /usr/bin/python3.13; do
    command -v "$p" >/dev/null 2>&1 && echo "  $p -> $("$p" --version 2>&1)"
done
echo "pip3=$(command -v pip3 2>&1)  uv=$(command -v uv 2>&1)"
echo "--- [debug] license dir BEFORE install ---"; ls -l "$MLROOT/licenses/" 2>&1 | sed 's/^/  /'

SD="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd)"
[ -f "$SD/prepare_run.m" ] || SD="/home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2026"
NNV_ROOT="$(cd "$SD/../../.." 2>/dev/null && pwd)"
[ -f "$NNV_ROOT/install.m" ] || NNV_ROOT="/home/ubuntu/toolkit/code/nnv"
REQ="$NNV_ROOT/tools/onnx2nnv_python/requirements.txt"
echo "SD=$SD  NNV_ROOT=$NNV_ROOT  REQ=$REQ ($([ -f "$REQ" ] && echo present || echo MISSING))"

# ---- 1) LICENSE (first: the MATLAB install in step 3 needs a licensed MATLAB) ----
echo "=== [1] LICENSE ==="
mkdir -p ~/.matlab/R2026a_licenses; cd ~/.matlab/R2026a_licenses
# MAC-locked R2026a license (HOSTID 02e1e896fadb == ENI eni-0b11771dfe21b94ee; exp 30-may-2027). URL QUOTED;
# bounded retry/time so a network problem fails fast instead of hanging (was --retry 100).
curl --retry 5 --retry-delay 3 --max-time 180 -fL -o license.lic "https://www.dropbox.com/scl/fi/w5jgddmf3qm5znjw67ajm/matlab-license-vnncomp2026-nnv.lic?rlkey=z3wnimbad4ykjyde95yhq7ik3&st=2oc64st3&dl=1"
echo "curl rc=$?"; sleep 2; ls -al license.lic 2>&1 | sed 's/^/  /'; echo "license head:"; head -3 license.lic 2>&1 | sed 's/^/  /'
if [ ! -s license.lic ] || ! grep -qE 'INCREMENT|MathWorks license|SERVER|DAEMON' license.lic; then echo "ERROR: license download failed/invalid"; fi
sudo mkdir -p "$MLROOT/licenses"
sudo cp -f license.lic "$MLROOT/licenses/" && echo "license copied to $MLROOT/licenses/" || echo "ERROR: license cp failed"
sudo rm -f "$MLROOT/licenses/license_info.xml"
echo "license dir AFTER copy:"; sudo ls -l "$MLROOT/licenses/" 2>&1 | sed 's/^/  /'

# ---- 2) system-python matlab.engine + ONNX stack + torch (sudo --break-system-packages -> system site-packages) ----
echo "=== [2] PYTHON: matlab.engine + ONNX stack into a SYSTEM python ==="
sudo apt-get install -y python3-pip software-properties-common 2>&1 | tail -4
ENGSRC="$MLROOT/extern/engines/python"
echo "engine source: $ENGSRC ($([ -d "$ENGSRC" ] && echo present || echo MISSING))"; ls -l "$ENGSRC" 2>&1 | sed 's/^/  /'
REQ_OK() { ( cd /tmp && "$1" -c "import matlab.engine,numpy,scipy,onnx,onnxruntime,onnxsim,onnxoptimizer,vnnlib" ) >/dev/null 2>&1; }
# Install engine FIRST; only if it imports do we spend time on the onnx stack + torch.
install_stack() {
    local PY="$1"
    echo "--- engine -> $PY ($("$PY" --version 2>&1)) ---"
    ( cd "$ENGSRC" && sudo "$PY" -m pip install . --break-system-packages ) 2>&1 | tail -20
    if ! ( cd /tmp && "$PY" -c "import matlab.engine" ) 2>/dev/null; then echo "  matlab.engine import FAILED on $PY -> skipping onnx/torch"; return 1; fi
    echo "  matlab.engine imports on $PY; installing onnx stack + torch"
    sudo "$PY" -m pip install -r "$REQ" --break-system-packages 2>&1 | tail -5
    sudo "$PY" -m pip install torch --break-system-packages 2>&1 | tail -3
    REQ_OK "$PY"
}
WORKING_PY=""; SEEN=""
for PY in /usr/bin/python3.13 /usr/bin/python3.12 /usr/bin/python3.11 /usr/bin/python3.10 /usr/bin/python3; do
    command -v "$PY" >/dev/null 2>&1 || continue
    RP="$(readlink -f "$PY")"; case " $SEEN " in *" $RP "*) continue;; esac; SEEN="$SEEN $RP"
    if install_stack "$PY"; then WORKING_PY="$PY"; echo "ENGINE_WORKS $PY"; break; fi
done
if [ -z "$WORKING_PY" ]; then
    echo "=== no pre-installed system python supports the R2026a engine (py3.10 too old); using deadsnakes ==="
    sudo add-apt-repository -y ppa:deadsnakes/ppa 2>&1 | tail -2
    sudo apt-get update 2>&1 | tail -2
    # CRITICAL (smoke 280): do NOT request python3.12-distutils -- distutils was removed in 3.12 so that
    # package does not exist, and naming it makes `apt-get install` abort the WHOLE transaction (python3.12
    # then never installs -> WORKING_SYS_PY=NONE). ensurepip is bundled in the stdlib for >=3.12. Try 3.12
    # then 3.11 then 3.13 (all R2026a-supported) and stop at the first whose matlab.engine imports.
    for V in 3.12 3.11 3.13; do
        sudo apt-get install -y "python$V" "python$V-venv" "python$V-dev" 2>&1 | tail -3
        PY="/usr/bin/python$V"
        command -v "$PY" >/dev/null 2>&1 || { echo "python$V not installed -> next"; continue; }
        echo "deadsnakes $PY -> $("$PY" --version 2>&1)"
        sudo "$PY" -m ensurepip --upgrade 2>&1 | tail -2
        if install_stack "$PY"; then WORKING_PY="$PY"; echo "ENGINE_WORKS $PY (deadsnakes $V)"; break; fi
    done
fi
# Seed the platform-prep python3 (prepare_instance runs onnx2nnv.py under bare python3).
sudo python3 -m pip install -r "$REQ" --break-system-packages 2>/dev/null || python3 -m pip install -r "$REQ" 2>/dev/null || true
echo "WORKING_SYS_PY=${WORKING_PY:-NONE}"

# ---- 3) NNV install + persist path -- NOW LICENSED. ONE `matlab -batch` (auto-quits on completion: no
#         lingering interactive MATLAB, no double-install, no path races). install.m does
#         tbxmanager(mpt/glpk/sedumi) + startup_nnv + self-test; then addpath+savepath persists NNV so the
#         per-instance runs (each a fresh matlab.engine session) resolve NNV via the saved path. We deliberately
#         DO NOT call prepare_run.m: it re-ran install (redundant ~2-3 min) and `disp(pyenv)` (MATLAB-calls-
#         python, which the engine-driven runs never use) was hanging toward the timeout (smoke 280). ----
echo "=== [3] NNV install + path persist (single matlab -batch, auto-quits) ==="
timeout 1200 matlab -batch "cd('${NNV_ROOT}'); install; addpath(genpath('${NNV_ROOT}')); savepath; fprintf('NNV_PATH_SAVED matlabroot=%s\n', matlabroot)" 2>&1 | tail -30
echo "NNV install+path rc=${PIPESTATUS[0]} (124=timeout)"
sudo chmod -R a+rX "$NNV_ROOT/tbxmanager" "$NNV_ROOT/code" 2>/dev/null || true

# ---- 4) GPU persistence + lock the driver ----
echo "=== [4] GPU ==="; nvidia-smi 2>&1 | head -5
sudo nvidia-smi -pm 1 2>/dev/null || true
NVPKG="$(dpkg -l 2>/dev/null | awk '/^ii +nvidia-driver-[0-9]/{print $2}' | head -1)"
sudo apt-mark hold linux-image-generic linux-headers-generic "${NVPKG:-nvidia-driver-570}" 2>/dev/null || true
sudo systemctl disable unattended-upgrades 2>/dev/null || true
echo "=== POST_INSTALL END $(date -u +%FT%TZ) WORKING_SYS_PY=${WORKING_PY:-NONE} MLROOT=$MLROOT ==="
