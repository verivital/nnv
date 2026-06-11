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

# MathWorks Package Manager: ONNX + PyTorch converters (needed by importNetworkFromONNX
# and the manifest importer). Add others here if the toolbox set grows.
apt-get update -y && apt-get install -y wget
wget -q https://www.mathworks.com/mpm/glnxa64/mpm
chmod +x mpm
./mpm install --release="${MATLAB_RELEASE}" \
    --products Deep_Learning_Toolbox_Converter_for_ONNX_Model_Format \
               Deep_Learning_Toolbox_Converter_for_PyTorch_Model_Format \
    --destination /usr/local/matlab || true

# One-time NNV toolbox install (tbxmanager: mpt/glpk/sedumi/...). Doing it here (not
# per-instance) keeps prepare/run overhead minimal. Adjust TOOLKIT if the image differs.
TOOLKIT="$(cd "$(dirname "$0")/../../../.." && pwd)"   # .../code/nnv
matlab -batch "cd('${TOOLKIT}'); install" || \
    echo "WARN: NNV install step failed here; run_instance will run startup_nnv as a fallback."

# Optional: Gurobi for a faster LP backend (uncomment + set license).
# cd ~ && wget -q https://packages.gurobi.com/12.0/gurobi12.0.0_linux64.tar.gz && tar xfz gurobi*.tar.gz
echo "Install complete."
