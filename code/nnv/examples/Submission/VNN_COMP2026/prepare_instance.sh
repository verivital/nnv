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
NNV_ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"   # VNN_COMP2026 -> .../code/nnv
gen_manifest() {
    local onnx="$1"
    local manifest="${onnx%.onnx}.nnv.mat"
    if [ ! -f "${manifest}" ]; then
        echo "Generating NNV manifest: ${manifest}"
        python3 "${NNV_ROOT}/tools/onnx2nnv_python/onnx2nnv.py" "${onnx}" "${manifest}" \
            || echo "WARN: manifest generation failed; run_instance will report error/unknown."
    fi
}
case "$2" in
    *lsnc_relu*|*traffic_signs*|*cgan*|*soundnessbench*|*vit*|*test*)
        gen_manifest "$3"
        ;;
    *nn4sys*)
        # Only the mscn models route via the manifest (lindex/pensieve use MATLAB's
        # importNetworkFromONNX). mscn_2048d_dual.onnx is corrupt upstream and abstains
        # to unknown in run_vnncomp_instance, so don't waste time trying to convert it.
        case "$3" in
            *mscn_2048d_dual*) : ;;
            *mscn*) gen_manifest "$3" ;;
        esac
        ;;
esac

exit 0
