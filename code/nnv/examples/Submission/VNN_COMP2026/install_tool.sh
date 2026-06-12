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
MATLAB_BIN="$(command -v matlab || echo /usr/local/matlab/bin/matlab)"
MATLAB_ROOT_RESOLVED="$(dirname "$(dirname "$(readlink -f "$MATLAB_BIN")")")"
pip3 install --no-cache-dir "${MATLAB_ROOT_RESOLVED}/extern/engines/python" || \
    pip3 install --no-cache-dir matlabengine || \
    echo "WARN: matlab.engine install failed; run_instance.sh will not work"
python3 -c "import matlab.engine" || echo "WARN: matlab.engine import check failed"

# Optional: Gurobi for a faster LP backend (uncomment + set license).
# cd ~ && wget -q https://packages.gurobi.com/12.0/gurobi12.0.0_linux64.tar.gz && tar xfz gurobi*.tar.gz
echo "Install complete."
