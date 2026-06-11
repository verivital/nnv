#!/bin/bash
# VNN-COMP 2026 prepare_instance.sh for NNV. Must be FAST (<=10 min budget) and do
# NO verification work -- only clean state so the per-instance run starts fresh.
# Args: v1, category, onnx_path, vnnlib_path.

TOOL_NAME="nnv"
VERSION_STRING="v1"
if [ "$1" != "${VERSION_STRING}" ]; then
    echo "Expected first argument (version string) '${VERSION_STRING}', got '$1'"
    exit 1
fi
echo "Preparing ${TOOL_NAME} for category '$2' (onnx='$3', vnnlib='$4')"

# Clear any zombie MATLAB/python from a previous instance so startup is clean.
killall -q python3 2>/dev/null
pkill -f -q matlab 2>/dev/null

exit 0
