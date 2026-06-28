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

# GLOBAL competition soundness-policy flags (NNV_CONV_TRUST_FP32 / NNV_CONV_NO_STAR / NNV_CONV_FRONTIER /
# NNV_QUARANTINE_CPSTAR) live in ONE place -- vnncomp2026_env.sh -- so every harness sources the SAME values
# and cannot drift (the 2026-06-24 cifar regression was sweep_lambda.sh dropping these). Full rationale +
# soundness notes are in that file.
source "$(dirname "${BASH_SOURCE[0]}")/vnncomp2026_env.sh"

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

# CONV alpha+beta-CROWN BaB: per-node beta split-duals (NNV_BAB_BETA_ITERS) optimized jointly with
# alpha, warm-started from the amortized root alpha (NNV_AMORT_ALPHA), refined over the still-
# undecided frontier in chunks of NNV_CONV_BETA_FRONTIER. This is the throughput lever that closes
# the loose root margin (cifar100 idx_2132 root -0.47 -> certified, 70s) the fixed-slope screen alone
# cannot. SOUNDNESS: beta>=0 is a valid Lagrangian dual on the active/inactive split constraints
# (the two children PARTITION the parent); keep-best + elementwise max(screen, refined) guarantees
# the refined margin is a valid, no-worse lower bound (parity + beta-Monte-Carlo gates green). A
# wrong split can only deepen the tree, never mis-certify. SCOPED to the conv resnet categories
# (cifar100/tinyimagenet/vggnet) -- the only ones whose Star never certifies and whose root margin
# needs beta; NOT set globally so FC/acas/falsify-primary categories are byte-identical. After
# spec-reduction these resnets keep ~1 unproven spec, so the autodiff tape is tiny and a large
# beta-frontier (64) batches the per-node refine cheaply.
case "$CATEGORY" in
    cifar100_2024)   # beta scoped to cifar100 ONLY: tinyimagenet (64x64) + vggnet beta time out the 100s budget (grind, no certs); env-tunable above for per-net experiments
        export NNV_BAB_BETA_ITERS=${NNV_BAB_BETA_ITERS:-3}
        export NNV_CONV_BETA_FRONTIER=${NNV_CONV_BETA_FRONTIER:-64}
        export NNV_AMORT_ALPHA=${NNV_AMORT_ALPHA:-20}
        ;;
esac

# Surface the post_install log here -- the platform's ToolkitPostInstall step log captures only its wrapper,
# so this is the only way to SEE why the python/venv setup did (or didn't) work on the eval box.
if [ -f "$HOME/.nnv_post_install.log" ]; then
    echo "=== POST_INSTALL LOG (begin) ==="; sed 's/^/[PI] /' "$HOME/.nnv_post_install.log"; echo "=== POST_INSTALL LOG (end) ==="
fi
echo "Running ${TOOL_NAME} on '$CATEGORY' (onnx='$ONNX_FILE', vnnlib='$VNNLIB_FILE', timeout=${TIMEOUT}s)"
# Run execute.py under the python that actually has matlab.engine (NNV_ORT_PYTHON, set by vnncomp2026_env.sh
# sourced above). A bare `python3` is anaconda on the eval box and lacks matlab.engine (smoke 269). GUARD: if
# the recorded python is missing/not executable, fall back to python3 (avoid the smoke-272 exit-127 dead end).
[ -x "$NNV_ORT_PYTHON" ] || NNV_ORT_PYTHON=python3
"${NNV_ORT_PYTHON:-python3}" "$EXECUTE" 'run_instance' "$CATEGORY" "$ONNX_FILE" "$VNNLIB_FILE" "$TIMEOUT" "$RESULTS_FILE"
