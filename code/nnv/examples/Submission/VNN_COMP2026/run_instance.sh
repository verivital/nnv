#!/bin/bash
# VNN-COMP 2026 run_instance.sh for NNV. Runs ONE instance and writes a single-word
# result (sat/unsat/unknown/timeout/error) + a validated SAT witness to RESULTS_FILE.
# Args: v1, category, onnx_path, vnnlib_path, results_file, timeout_seconds.
#
# Self-contained 2026 submission: uses the LOCAL run_vnncomp_instance.m (gradient PGD
# falsifier + witness validation) via the local execute.py bridge. The 2026 folder owns
# its full verification logic; VNN_COMP2025 is frozen at its original 2025 submission.

TOOL_NAME="nnv"
VERSION_STRING="v1"
if [ "$1" != "${VERSION_STRING}" ]; then
    echo "Expected first argument (version string) '${VERSION_STRING}', got '$1'"
    exit 1
fi
CATEGORY="$2"; ONNX_FILE="$3"; VNNLIB_FILE="$4"; RESULTS_FILE="$5"; TIMEOUT="$6"

DIR="$(cd "$(dirname "$0")" && pwd)"
EXECUTE="$DIR/execute.py"   # local maintained MATLAB<->result bridge (self-contained 2026)
if [ ! -f "$EXECUTE" ]; then
    echo "error" > "$RESULTS_FILE"
    echo "ERROR: execute.py not found at $EXECUTE"; exit 1
fi

# Conv GPU-BaB tuning (FIX A1/A3) -- exported into the MATLAB process via execute.py.
# A1: default the slow UNSOUND FP32 GPU screen OFF so the cheap SOUND FP64-CPU confirm runs DIRECTLY
#     (on hard cifar/tinyimagenet/challenging instances the screen never returns 'robust' in budget,
#     starving the confirm that certifies in <31 nodes). The emit still comes from the FP64 double
#     pass + argmax-spec + orientation guards -> sound; this only changes WHICH sound path runs first.
# A3: widen the conv BaB frontier to fight the launch-bound tail (default 32). A frontier/node cap can
#     only yield 'unknown', never a wrong verdict. Both vars affect ONLY the conv precheck path;
#     non-conv categories ignore them. (NNV_SOUND_FP32_TIGHT / FIX A2 intentionally NOT set yet -- it
#     needs the owed known-SAT conv soundness test first.)
export NNV_CONV_GPU_SCREEN=0
export NNV_CONV_FRONTIER=256

echo "Running ${TOOL_NAME} on '$CATEGORY' (onnx='$ONNX_FILE', vnnlib='$VNNLIB_FILE', timeout=${TIMEOUT}s)"
python3 "$EXECUTE" 'run_instance' "$CATEGORY" "$ONNX_FILE" "$VNNLIB_FILE" "$TIMEOUT" "$RESULTS_FILE"
