#!/bin/bash
# VNN-COMP 2026 run_instance.sh for NNV. Runs ONE instance and writes a single-word
# result (sat/unsat/unknown/timeout/error) + a validated SAT witness to RESULTS_FILE.
# Args: v1, category, onnx_path, vnnlib_path, results_file, timeout_seconds.
#
# Reuses the maintained MATLAB runner in ../VNN_COMP2025 (run_vnncomp_instance.m, now
# with the gradient PGD falsifier + witness validation) via the existing execute.py
# bridge, so there is a single source of truth for the verification logic.

TOOL_NAME="nnv"
VERSION_STRING="v1"
if [ "$1" != "${VERSION_STRING}" ]; then
    echo "Expected first argument (version string) '${VERSION_STRING}', got '$1'"
    exit 1
fi
CATEGORY="$2"; ONNX_FILE="$3"; VNNLIB_FILE="$4"; RESULTS_FILE="$5"; TIMEOUT="$6"

DIR="$(cd "$(dirname "$0")" && pwd)"
EXECUTE="$DIR/../VNN_COMP2025/execute.py"   # maintained MATLAB<->result bridge
if [ ! -f "$EXECUTE" ]; then
    echo "error" > "$RESULTS_FILE"
    echo "ERROR: execute.py not found at $EXECUTE"; exit 1
fi

echo "Running ${TOOL_NAME} on '$CATEGORY' (onnx='$ONNX_FILE', vnnlib='$VNNLIB_FILE', timeout=${TIMEOUT}s)"
python3 "$EXECUTE" 'run_instance' "$CATEGORY" "$ONNX_FILE" "$VNNLIB_FILE" "$TIMEOUT" "$RESULTS_FILE"
