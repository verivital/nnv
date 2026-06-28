#!/bin/bash
# VNN-COMP 2026 install_tool.sh for NNV (https://github.com/verivital/nnv).
# Runs ONCE on the evaluation image (Ubuntu, MATLAB R2026a), as the NON-ROOT toolkit user
# (config.yaml run_installation_script_as_root: False).
#
# MINIMAL BY DESIGN (matches the known-good 2025 pattern). This script ONLY installs the MATLAB
# support packages via mpm, AS THE RUNNING USER -- so they land in ~/Documents/MATLAB/SupportPackages
# and stay visible to the toolkit user that runs verification (run_toolkit_as_root: False). If this
# ran as root, the support packages would relocate to /root/... and the ONNX converter would vanish.
#
# Everything that needs root (apt) or a LICENSED MATLAB (the NNV tbxmanager install) lives in
# post_install.sh, which runs AFTER licensing. The 2026-06-28 smoke test (task 266) proved why:
# doing `apt-get install` here failed "Permission denied (are you root?)" (non-root), and
# `matlab -batch install` failed license Error -1.2 (the license isn't installed until post_install).
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
echo "Installing ${TOOL_NAME} support packages (MATLAB ${MATLAB_RELEASE})"

ip link show   # mac address (licensing)
echo "$USER"    # username (licensing)
mkdir -p "$HOME/.matlab/${MATLAB_RELEASE}_licenses"

# ---- detect the MATLAB install root (DO NOT hardcode): works for the MathWorks AMI
# (/usr/local/matlab) AND a self-installed user-dir install. Explicit MATLABROOT env wins. ----
if [ -n "${MATLABROOT:-}" ]; then
    :
elif MATLAB_BIN_EARLY="$(command -v matlab)"; then
    MATLAB_REAL_EARLY="$(readlink -f "$MATLAB_BIN_EARLY" 2>/dev/null || echo "$MATLAB_BIN_EARLY")"
    MATLABROOT="$(dirname "$(dirname "$MATLAB_REAL_EARLY")")"
else
    MATLABROOT="/usr/local/matlab"
fi
if [ ! -x "${MATLABROOT}/bin/matlab" ]; then
    echo "WARN: MATLABROOT='${MATLABROOT}' has no bin/matlab; set MATLABROOT explicitly if installs land wrong."
fi
echo "Using MATLABROOT=${MATLABROOT}"

# MathWorks Package Manager (mpm): provision the toolbox set NNV needs PLUS the ONNX/PyTorch converters,
# AS THE USER (default destination ~/Documents/MATLAB/SupportPackages so they're visible to the toolkit
# user). mpm is ADDITIVE + idempotent. Override via env NNV_MPM_PRODUCTS.
# wget is needed to fetch mpm; it is usually preinstalled on the MathWorks AMI, so the apt is best-effort
# (sudo works on this platform per the 2025 post_install, but guard it so a non-root/no-net box can't abort).
sudo apt-get install -y wget 2>/dev/null || apt-get install -y wget 2>/dev/null || echo "WARN: could not apt-install wget (assuming preinstalled)"
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

# ---- GPU driver / CUDA sanity (warn only; post_install holds the driver). R2026a (CUDA 12.8) needs >=570. ----
if command -v nvidia-smi >/dev/null 2>&1; then
    DRV="$(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1)"
    DRV_MAJ="${DRV%%.*}"
    if [ -z "$DRV" ] || ! [ "$DRV_MAJ" -eq "$DRV_MAJ" ] 2>/dev/null; then
        echo "WARN: nvidia-smi present but driver version undetermined; GPU compute likely unavailable (need >=570)."
    elif [ "$DRV_MAJ" -lt 570 ]; then
        echo "WARN: NVIDIA driver ${DRV} < 570 -> R2026a (CUDA 12.8) GPU compute will FAIL."
    else
        echo "NVIDIA driver ${DRV} OK for R2026a CUDA 12.8 (>=570)."
    fi
else
    echo "NOTE: nvidia-smi not found -> GPU-BaB CPU-only. Install NVIDIA driver >=570 for GPU."
fi

echo "Support-package install complete. License + licensed MATLAB setup + apt/pip deps run in post_install.sh."
