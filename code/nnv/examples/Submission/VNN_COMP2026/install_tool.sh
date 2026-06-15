#!/bin/bash
# VNN-COMP 2026 install_tool.sh for NNV (https://github.com/verivital/nnv).
# Runs ONCE on the evaluation image (Ubuntu, MATLAB R2026a). Installs the MATLAB
# support packages NNV needs and warms the tbxmanager toolboxes so per-instance
# overhead stays low (NNV had the worst startup overhead in 2025 -- 14.2s).
#
# Arg 1: version string (must match VERSION_STRING below).

set -e
TOOL_NAME="nnv"
VERSION_STRING="v1"
MATLAB_RELEASE="R2026a"

if [ "$1" != "${VERSION_STRING}" ]; then
    echo "Expected first argument (version string) '${VERSION_STRING}', got '$1'"
    exit 1
fi
echo "Installing ${TOOL_NAME} dependencies (MATLAB ${MATLAB_RELEASE})"

ip link show   # mac address (licensing)
echo "$USER"    # username (licensing)
mkdir -p "$HOME/.matlab/${MATLAB_RELEASE}_licenses"

# ---- detect the MATLAB install root (DO NOT hardcode): works for the MathWorks AMI
# (/usr/local/matlab) AND a self-installed user-dir install (e.g. ~/MATLAB/R2026a). An
# explicit MATLABROOT env wins; otherwise derive from `matlab` on PATH; else fall back. ----
if [ -n "${MATLABROOT:-}" ]; then
    :   # explicit override -- trust the caller
elif MATLAB_BIN_EARLY="$(command -v matlab)"; then
    MATLAB_REAL_EARLY="$(readlink -f "$MATLAB_BIN_EARLY" 2>/dev/null || echo "$MATLAB_BIN_EARLY")"
    MATLABROOT="$(dirname "$(dirname "$MATLAB_REAL_EARLY")")"
else
    MATLABROOT="/usr/local/matlab"
fi
# Sanity: a real MATLAB root has bin/matlab. If not (e.g. `matlab` was a wrapper in /usr/local/bin
# so dirname^2 mis-resolved), warn -- mpm/pip would otherwise target the wrong place.
if [ ! -x "${MATLABROOT}/bin/matlab" ]; then
    echo "WARN: MATLABROOT='${MATLABROOT}' has no bin/matlab; set MATLABROOT explicitly if installs land wrong."
fi
echo "Using MATLABROOT=${MATLABROOT}"

# MathWorks Package Manager (mpm): provision the FULL MATLAB toolbox set NNV needs PLUS the
# ONNX/PyTorch converters, into the detected MATLAB root. mpm is ADDITIVE + idempotent --
# products already present are skipped -- so this is fast on the AMI (base toolboxes
# preinstalled) AND completes a general/from-scratch install. Override via env NNV_MPM_PRODUCTS.
apt-get update -y && apt-get install -y wget
wget -q https://www.mathworks.com/mpm/glnxa64/mpm
chmod +x mpm
NNV_MPM_PRODUCTS="${NNV_MPM_PRODUCTS:-Deep_Learning_Toolbox Parallel_Computing_Toolbox \
Computer_Vision_Toolbox Statistics_and_Machine_Learning_Toolbox Optimization_Toolbox \
Control_System_Toolbox Image_Processing_Toolbox Symbolic_Math_Toolbox \
System_Identification_Toolbox Deep_Learning_Toolbox_Converter_for_ONNX_Model_Format \
Deep_Learning_Toolbox_Converter_for_PyTorch_Model_Format}"
./mpm install --release="${MATLAB_RELEASE}" \
    --destination "${MATLABROOT}" \
    --accept-vendor-licenses \
    --products ${NNV_MPM_PRODUCTS} || \
    echo "WARN: mpm install returned non-zero (some products may already be present)"

# One-time NNV toolbox install (tbxmanager: mpt/glpk/sedumi/...). Doing it here (not
# per-instance) keeps prepare/run overhead minimal. Adjust TOOLKIT if the image differs.
TOOLKIT="$(cd "$(dirname "$0")/../../.." && pwd)"   # VNN_COMP2026 -> .../code/nnv (has install.m)
matlab -batch "cd('${TOOLKIT}'); install" || \
    echo "WARN: NNV install step failed here; run_instance will run startup_nnv as a fallback."

# Python importer + onnxruntime. onnx2nnv.py (tools/onnx2nnv_python) converts the models
# MATLAB's importer cannot parse (lsnc_relu, traffic_signs, cgan, soundnessbench) into NNV
# manifests in prepare_instance.sh; onnxruntime also backs the SAT-witness replay gate that
# prevents false-SAT (-150) verdicts. Pin a known-good set.
apt-get install -y python3 python3-pip
pip3 install --no-cache-dir onnx==1.20.0 onnxruntime==1.23.1 numpy scipy
# onnx2nnv.py (the manifest importer for lsnc_relu/traffic_signs/cgan/soundnessbench)
# hard-requires the simplifier stack; without these prepare_instance.sh fails to
# generate manifests and every manifest-category instance errors (2026-06-12 dry run).
pip3 install --no-cache-dir onnxsim onnxoptimizer

# MATLAB Engine API for Python: execute.py (the run_instance.sh bridge) does
# `import matlab.engine`; without it EVERY instance fails instantly (caught in
# the 2026-06-12 AWS dry run on the MathWorks R2026a AMI). Install from the
# MATLAB tree -- guaranteed version match with the installed release; the PyPI
# matlabengine wheel can fail to build and may mismatch the MATLAB version.
pip3 install --no-cache-dir "${MATLABROOT}/extern/engines/python" || \
    pip3 install --no-cache-dir matlabengine || \
    echo "WARN: matlab.engine install failed (import check below decides)"
# Fatal gate: execute.py does `import matlab.engine` on every instance, so a
# missing engine means zero points. Fail the install loudly rather than let a
# misleading "Install complete." defer the failure to runtime.
python3 -c "import matlab.engine" || { echo "ERROR: matlab.engine import failed; fix before running instances" >&2; exit 1; }

# ---- GPU driver / CUDA sanity (NNV's GPU-BaB engine needs a working gpuDevice) ------------
# R2026a ships CUDA 12.8 -> requires NVIDIA driver >= 570. Warn (do NOT fail) when the GPU is
# unusable so the user can upgrade (sudo apt install nvidia-driver-580 && sudo reboot); CPU
# verification still works without a GPU.
if command -v nvidia-smi >/dev/null 2>&1; then
    DRV="$(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1)"
    DRV_MAJ="${DRV%%.*}"
    if [ -z "$DRV" ] || ! [ "$DRV_MAJ" -eq "$DRV_MAJ" ] 2>/dev/null; then
        echo "WARN: nvidia-smi present but the driver version could not be determined (driver not loaded / no GPU?)."
        echo "      GPU compute (gpuDevice) likely unavailable; R2026a (CUDA 12.8) needs NVIDIA driver >=570."
    elif [ "$DRV_MAJ" -lt 570 ]; then
        echo "WARN: NVIDIA driver ${DRV} < 570 -> R2026a (CUDA 12.8) GPU compute (gpuDevice) will FAIL."
        echo "      Fix: sudo apt install -y nvidia-driver-580 && sudo reboot"
    else
        echo "NVIDIA driver ${DRV} OK for R2026a CUDA 12.8 (requires >=570)."
    fi
else
    echo "NOTE: nvidia-smi not found -> GPU-BaB will run CPU-only. Install NVIDIA driver >=570 for GPU."
fi

# Optional: Gurobi for a faster LP backend (uncomment + set license).
# cd ~ && wget -q https://packages.gurobi.com/12.0/gurobi12.0.0_linux64.tar.gz && tar xfz gurobi*.tar.gz
echo "Install complete."
