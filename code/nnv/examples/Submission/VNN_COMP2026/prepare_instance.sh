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

# Manifest categories: MATLAB's importNetworkFromONNX cannot parse these models, so NNV
# imports them via the Python importer (onnx2nnv.py -> <model>.nnv.mat). Generating the
# manifest here is allowed setup -- it is format CONVERSION, not analysis -- and is cheap
# (~1-3 s/model, far under the 10 min prepare cap). The forward pass of each manifest is
# cross-validated vs onnxruntime; any SAT witness is replayed through onnxruntime before
# being emitted, so a mis-import degrades to error/unknown, never an unsound verdict.
case "$2" in
    *lsnc_relu*|*traffic_signs*|*cgan*|*soundnessbench*|*vit*)
        NNV_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"   # VNN_COMP2026 -> .../code/nnv
        MANIFEST="${3%.onnx}.nnv.mat"
        if [ ! -f "${MANIFEST}" ]; then
            echo "Generating NNV manifest: ${MANIFEST}"
            python3 "${NNV_ROOT}/tools/onnx2nnv_python/onnx2nnv.py" "$3" "${MANIFEST}" \
                || echo "WARN: manifest generation failed; run_instance will report error/unknown."
        fi
        ;;
esac

exit 0
