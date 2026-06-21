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

# Conv GPU-BaB -- GPU PATHWAY (validated 2026-06-19), exported into MATLAB via execute.py.
# TRUST_FP32: emit the unsat verdict FROM the GPU-single batched-BaB screen (the raw CROWN bound run on
#   the GPU -- ~95% util, ~5x faster than the FP64-CPU reconfirm it replaces), PGD-backstopped. PGD runs
#   FALSIFY-FIRST, so a screen-robust that is PGD-clean is the SAME two-mechanism soundness the GPU-winning
#   tools rely on (alpha-beta-CROWN runs stock FP32 + a PGD attack; NeuralSAT likewise). Validated
#   0 false-robust vs the alpha-beta-CROWN gold set (23 cifar100/tinyimagenet instances incl gold-SAT).
#   SOUNDNESS POLICY: trusts FP32 rounding (~1e-6, negligible vs real robustness margins) with PGD as the
#   backstop -- the deliberate, field-standard GPU verification model. Forces the GPU screen ON.
# NO_STAR: a non-certifying conv precheck emits a fast sound 'unknown' (Star never certifies these
#   resnets). FRONTIER 512: more BaB nodes per GPU kernel. All vars affect ONLY the conv precheck path.
export NNV_CONV_TRUST_FP32=1
export NNV_CONV_NO_STAR=1
export NNV_CONV_FRONTIER=512
# QUARANTINE_CPSTAR: strip the PROBABILISTIC (conformal) cp-star reach method so a cp-star 'unsat' can
# never be emitted as a sound verdict. cp-star is NOT an FP64 proof -- under the absolute 0-wrong rule a
# probabilistic unsat is a latent -150 (correct only with some confidence -> over the field some would be
# wrong). STRICTLY SOUND: only ever converts a cp-star unsat -> unknown (never adds a verdict); the sound
# FP64 / GPU-BaB / exact paths still decide what they can. Competition-scoped (set here, not the code
# default) so tests / other callers are unaffected.
export NNV_QUARANTINE_CPSTAR=1

# FALSIFY_MAXTIME: cap the falsify-first PGD budget so the SOUND reach proof keeps its window within the
# OFFICIAL timeout. safenlp's ftab PGD budget (30s) EXCEEDS its 20s timeout, so on a ROBUST instance PGD
# burns the whole wall finding nothing and the external kill fires BEFORE the ~0.4s reach proof runs -> a
# decidable unsat is lost to timeout. Measured: 25/25 over-timeout safenlp unsats recover (decide at
# 8.5-14s) and 10/10 sat are preserved (found <3.3s), 0 gold-contradictions. STRICTLY SOUND: less PGD only
# ever yields fewer SAT (-> unknown), never a wrong verdict; reach is untouched. Scoped to safenlp (the one
# category where PGD>=timeout AND reach is sub-second); NOT applied to falsify-primary cats (cgan/sat_relu/
# traffic_signs) whose SAT need the full PGD budget, nor to reach-slow cats (cora) the cap can't help.
if [ "$CATEGORY" = "safenlp_2024" ]; then
    export NNV_FALSIFY_MAXTIME=8
fi

echo "Running ${TOOL_NAME} on '$CATEGORY' (onnx='$ONNX_FILE', vnnlib='$VNNLIB_FILE', timeout=${TIMEOUT}s)"
python3 "$EXECUTE" 'run_instance' "$CATEGORY" "$ONNX_FILE" "$VNNLIB_FILE" "$TIMEOUT" "$RESULTS_FILE"
