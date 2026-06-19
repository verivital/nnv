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
    *lsnc_relu*|*traffic_signs*|*cgan*|*soundnessbench*|*vit*|*test*|*linearize*)
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

# Net cache: importNetworkFromONNX + matlab2nnv is ONNX-determined, so pre-build the NNV net cache
# ONCE per onnx HERE (untimed prepare; VNN-COMP allows ONNX->format conversion in prepare, same as
# the manifests above) so the timed run loads a cached net (.netcache.mat) instead of re-converting
# (~10 s/instance saved; a category with many vnnlib instances pays the conversion once). Skipped
# if the cache already exists. Soundness: the cache is keyed (onnx size+mtime+NNV version) so a hit
# is only for the SAME onnx -- run_vnncomp_instance re-imports on any mismatch. Pre-build via
# NNV_PREP_CACHE=1, which makes run_vnncomp_instance return right after building the cache.
MATLAB_BIN="${MATLAB_BIN:-matlab}"
build_netcache() {
    local cat="$1" onnx="$2" vnnlib="$3"
    # Name MUST match the MATLAB side exactly: i_load_vnncomp_network_cached writes/reads
    # [char(onnx) '.netcache.mat'] -> "<onnx>.netcache.mat" (APPENDS; does not strip .onnx).
    # A stripped name ("model.netcache.mat") would never match, so this skip-guard would never
    # fire and prepare would relaunch MATLAB on every instance instead of once per onnx.
    local cache="${onnx}.netcache.mat"
    [ -f "${cache}" ] && return
    echo "Pre-building NNV net cache: ${cache}"
    NNV_PREP_CACHE=1 "${MATLAB_BIN}" -batch \
        "cd('${NNV_ROOT}'); startup_nnv; addpath('${NNV_ROOT}/examples/Submission/VNN_COMP2026'); cd('${NNV_ROOT}/examples/Submission/VNN_COMP2026'); try, run_vnncomp_instance('${cat}','${onnx}','${vnnlib}',tempname); catch ME, fprintf(2,'%s\n',ME.message); end" \
        >/dev/null 2>&1 || echo "WARN: net cache pre-build failed; run_instance converts on first use."
}
# matlab2nnv (non-manifest) categories whose single (or few) network(s) are reused across many
# vnnlib instances. Manifest categories above already cache via the .nnv.mat.
case "$2" in
    *cifar100*|*tinyimagenet*|*vggnet*|*yolo*|*metaroom*|*malbeware*|*collins_rul*|*dist_shift*|*relusplitter*|*safenlp*|*sat_relu*|*cersyve*|*tllverify*|*cora*|*acasxu*|*challenging*)
        build_netcache "$2" "$3" "$4" ;;
esac

exit 0
