#!/usr/bin/env bash
# Run the full NNV3.0 ATVA 2026 artifact-evaluation suite (FairNNV +
# ProbVer + GNNV + VideoStar + ModelStar + ToolComparison) in one shot.
# Each experiment runs in its own MATLAB session so a crash in one does
# not lose the others.
#
# Usage (inside the nnv3.0 container, or any MATLAB+NNV install):
#   bash code/nnv/examples/NNV3.0/run_all.sh
#
# Optional environment overrides:
#   NNV3_SKIP="probver videostar"     # space-separated experiments to skip
#                                       (toolcomparison is also a valid name)
#   NNV3_LOG_DIR=/tmp/nnv3_logs       # where per-experiment stdout lands
#   MATLAB=/usr/local/bin/matlab      # MATLAB binary (default: 'matlab' on PATH)
#   TOOLCOMPARISON_MODE=smoke|full    # ToolComparison mode (default: smoke).
#                                       'full' takes ~3-5 h and renders the
#                                       ATVA 2026 paper's Tables 5, 6, 7.
#
# Exit code is 0 only if every selected experiment finished cleanly. The final
# summary table (and a CSV at $NNV3_LOG_DIR/summary.csv) lists wall-clock time
# per experiment.

set -uo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    cat <<'EOF'
Usage: bash run_all.sh

Runs the NNV3.0 repeatability suite end-to-end (FairNNV + ProbVer + GNNV
+ VideoStar + ModelStar + ToolComparison). Each experiment runs in its
own MATLAB session so a crash in one does not lose the others.

Environment overrides:
  NNV3_SKIP="probver videostar"   space-separated experiment names to skip
                                  (valid: fairnnv probver gnnv videostar
                                  modelstar toolcomparison)
  NNV3_LOG_DIR=/tmp/nnv3_logs     where per-experiment stdout lands
                                  (default: <script_dir>/repeatability_logs)
  MATLAB=/usr/local/bin/matlab    MATLAB binary (default: 'matlab' on PATH)
  NNV3_FORCE_GPU=0|1              override GPU autodetection
  NNV3_FORCE_MEMORY=1             attempt ProbVer regardless of host RAM
  NNV3_MIN_MEMORY_GB=48           ProbVer auto-skip threshold (default 48)
  TOOLCOMPARISON_MODE=smoke|full  ToolComparison mode (default: smoke,
                                  ~12 min). 'full' takes ~3-5 h and renders
                                  the ATVA 2026 paper's Tables 5, 6, and 7.

Exit code is 0 only if every selected experiment finished cleanly. A
summary CSV at $NNV3_LOG_DIR/summary.csv lists wall-clock time per
experiment.
EOF
    exit 0
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NNV_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"  # .../code/nnv
LOG_DIR="${NNV3_LOG_DIR:-${SCRIPT_DIR}/repeatability_logs}"
MATLAB="${MATLAB:-matlab}"
SKIP="${NNV3_SKIP:-}"

mkdir -p "$LOG_DIR"
SUMMARY_CSV="${LOG_DIR}/summary.csv"
echo "experiment,status,wall_seconds,log" > "$SUMMARY_CSV"

# Single consolidated log captures everything the terminal sees.
RUN_LOG="${LOG_DIR}/run.log"
: > "$RUN_LOG"

# Terminal filter: keeps only status markers, per-instance verdicts, milestone
# summaries, and final result tables. Everything else is dropped from both the
# terminal stream and run.log. Set NNV3_VERBOSE=1 to disable the filter (every
# line passes through unchanged, for debugging).
TERMINAL_WHITELIST='^=== |^\[run_all\]|^\[install|^\[make_table_main\]|^\[OK\]|^\[AIVL\]|^======= |^Model: |^Number of fair |^Number of non-fair |^Number of unknown |^It took a total of|idx=[0-9]+ status=|^Fraction:|Percentage of images verified|verified=[0-9]+/[0-9]+|^Starting verification with epsilon|^Epsilon = [0-9]+/255:|^  (Verified|Violated|Unknown|Timeout|Average time):|verified \([0-9.]+|violated \([0-9.]+|timeout \(>?[0-9.]+|unknown \([0-9.]+|^Trained parameters saved|^Final loss:|^Estimated Lipschitz|^Using device:|sufficient amount of principal directions|directions and training time saved|^Benchmark .*Tool .*Algorithm|^-{10,}$|^(acas_xu_p[34]|rl|oval21|collins_rul|mnist_resnet8) +(nnv|aivl)|^=== run_toolcomparison|^(experiment|fairnnv|probver|gnnv|videostar|modelstar|toolcomparison)[, ]|^Consolidated log:|^Results consolidated'

filter_terminal() {
    if [[ "${NNV3_VERBOSE:-0}" == "1" ]]; then
        cat
    else
        grep -E --line-buffered "$TERMINAL_WHITELIST" || true
    fi
}

# ---------------------------------------------------------------------------
# GPU autodetection. ProbVer's cp-star reachability trains a surrogate
# network on CUDA via the Python venv; without a GPU it cannot run, so we
# auto-skip it. FairNNV and ModelStar are CPU-only; GNNV and VideoStar
# use MATLAB's reachability which is largely CPU but with some
# GPU-accelerated paths — they will run on a CPU host but expect
# significant slowdowns. Override the autodetection with
# NNV3_FORCE_GPU={0,1}.
# ---------------------------------------------------------------------------
detect_gpu() {
    if [[ -n "${NNV3_FORCE_GPU:-}" ]]; then
        [[ "${NNV3_FORCE_GPU}" == "1" ]]
        return $?
    fi
    if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1 \
       && nvidia-smi -L | grep -q '^GPU '; then
        return 0
    fi
    return 1
}

if detect_gpu; then
    HAS_GPU=1
    echo "[run_all] GPU detected:"
    nvidia-smi -L 2>/dev/null | sed 's/^/[run_all]   /'
else
    HAS_GPU=0
    echo "[run_all] No NVIDIA GPU detected — deferring to CPU."
    echo "[run_all]   ProbVer will be SKIPPED (cp-star reachability requires CUDA)."
    echo "[run_all]   FairNNV is CPU-only, no change."
    echo "[run_all]   GNNV / VideoStar will attempt CPU-only execution;"
    echo "[run_all]   expect significant slowdown vs the GPU baselines in the README."
    # Mark probver as skipped via the same SKIP mechanism users have.
    if [[ " $SKIP " != *" probver "* ]]; then
        SKIP="${SKIP} probver"
    fi
fi
export NNV3_HAS_GPU="$HAS_GPU"

# ---------------------------------------------------------------------------
# Memory preset. ProbVer's cp-star reachability builds ImageStars whose peak
# size on TinyYOLO has been observed to exceed 31 GB on this host, causing a
# SIGKILL inside `Prob_reach` (see ProbReach_ImageStar warning "The Image
# Star is large for your memory and should be presented in sparse format").
# We require >= NNV3_MIN_MEMORY_GB (default 48) of container RAM for ProbVer
# and auto-skip it below that, unless NNV3_FORCE_MEMORY=1.
# FairNNV / GNNV / VideoStar / ModelStar fit comfortably in <= 16 GB.
# ---------------------------------------------------------------------------
MIN_MEMORY_GB="${NNV3_MIN_MEMORY_GB:-48}"
MEM_TOTAL_GB="$(awk '/MemTotal/ { printf "%d", $2/1024/1024 }' /proc/meminfo 2>/dev/null || echo 0)"
echo "[run_all] Container RAM: ${MEM_TOTAL_GB} GB (recommended >= ${MIN_MEMORY_GB} GB for full suite)."
if (( MEM_TOTAL_GB < MIN_MEMORY_GB )); then
    if [[ "${NNV3_FORCE_MEMORY:-0}" == "1" ]]; then
        echo "[run_all]   NNV3_FORCE_MEMORY=1 set — running ProbVer anyway; OOM risk on instance 52/68."
    elif [[ " $SKIP " == *" probver "* ]]; then
        : # already skipped
    else
        echo "[run_all]   ProbVer auto-SKIPPED (memory below threshold)."
        echo "[run_all]   On Windows hosts, raise the WSL2 cap by creating ~/.wslconfig with:"
        echo "[run_all]       [wsl2]"
        echo "[run_all]       memory=56GB"
        echo "[run_all]       swap=32GB"
        echo "[run_all]   then run 'wsl --shutdown' and restart Docker Desktop."
        echo "[run_all]   On Linux, add memory if available; otherwise set NNV3_FORCE_MEMORY=1 to try anyway."
        SKIP="${SKIP} probver"
        SKIP_REASON_PROBVER="memory ${MEM_TOTAL_GB}GB < ${MIN_MEMORY_GB}GB threshold"
    fi
fi

# MATLAB prelude executed at the start of every experiment's MATLAB session.
# Forward-compat call (no-op on hosts without an NVIDIA GPU) keeps Blackwell /
# RTX 50-series usable. Warning silencing strips four known-noisy warnings
# from the terminal/log stream:
#   - nnet_cnn_onnx:onnx:WarnAPIDeprecation
#       "'importONNXNetwork' is not recommended..." (R2025b deprecation of an
#       API the artifact deliberately uses; no actionable signal for reviewers).
#   - nnet_cnn_onnx:onnx:InputLabelMismatch
#       "Data format '...' you specified for input 1 does not match the format
#       '...' derived by the software." (Our ONNX nets pass an explicit input
#       format; MATLAB falls back to its derived format which is identical in
#       intent.)
#   - nnet_cnn:internal:cnn:analyzer:NetworkAnalyzer:NetworkHasWarnings
#       Post-import analyzer warnings about cosmetic layer-property defaults
#       (e.g., "Empty Classes property") that don't affect reachability.
#   - MATLAB:mpath:nameNonexistentOrNotADirectory
#       Defensive: addpath() on a non-existent dir. The .dockerignore now
#       excludes stale runtime-output dirs (repeatability_logs, etc.) so
#       pathdef.m no longer hard-codes refs to them, but this catches any
#       transient saved-path entry that survives a future build.
# `warning('off','backtrace')` drops the multi-line stack trace each remaining
# warning would otherwise emit (each trace adds ~8 lines per warning event).
MATLAB_PRELUDE="\
warning('off','backtrace'); \
warning('off','nnet_cnn_onnx:onnx:WarnAPIDeprecation'); \
warning('off','nnet_cnn_onnx:onnx:InputLabelMismatch'); \
warning('off','nnet_cnn:internal:cnn:analyzer:NetworkAnalyzer:NetworkHasWarnings'); \
warning('off','MATLAB:mpath:nameNonexistentOrNotADirectory'); \
addpath(genpath('${NNV_ROOT}')); \
try, parallel.gpu.enableCUDAForwardCompatibility(true); catch; end;"

run_one() {
    # Args: name, subdir, kind ("matlab"|"shell"), entry
    local name="$1"
    local subdir="$2"
    local kind="$3"
    local entry="$4"

    if [[ " $SKIP " == *" $name "* ]]; then
        local reason="NNV3_SKIP"
        if [[ "$name" == "probver" ]]; then
            if [[ "$HAS_GPU" -eq 0 ]]; then
                reason="no GPU (CPU fallback unsupported for cp-star)"
            elif [[ -n "${SKIP_REASON_PROBVER:-}" ]]; then
                reason="$SKIP_REASON_PROBVER"
            fi
        fi
        printf '\n=== %-12s SKIPPED (%s) ===\n' "$name" "$reason" | tee -a "$RUN_LOG"
        echo "${name},skipped,0,${reason}" >> "$SUMMARY_CSV"
        return 0
    fi

    printf '\n=== %-12s start: %s ===\n' "$name" "$(date -u +%FT%TZ)" | tee -a "$RUN_LOG"
    local t0
    t0=$(date -u +%s)
    if [[ "$kind" == "matlab" ]]; then
        (
            cd "${SCRIPT_DIR}/${subdir}" || exit 99
            "$MATLAB" -nodisplay -batch "${MATLAB_PRELUDE} run('${entry}'); exit()"
        ) 2>&1 | filter_terminal | tee -a "$RUN_LOG"
    else
        (
            cd "${SCRIPT_DIR}/${subdir}" || exit 99
            bash "${SCRIPT_DIR}/${subdir}/${entry}"
        ) 2>&1 | filter_terminal | tee -a "$RUN_LOG"
    fi
    local rc=${PIPESTATUS[0]}
    local t1
    t1=$(date -u +%s)
    local elapsed=$(( t1 - t0 ))
    if [[ $rc -eq 0 ]]; then
        printf '=== %-12s OK in %ds ===\n' "$name" "$elapsed" | tee -a "$RUN_LOG"
        echo "${name},ok,${elapsed},${RUN_LOG}" >> "$SUMMARY_CSV"
    else
        printf '=== %-12s FAILED (rc=%d) in %ds ===\n' "$name" "$rc" "$elapsed" | tee -a "$RUN_LOG"
        echo "${name},failed,${elapsed},${RUN_LOG}" >> "$SUMMARY_CSV"
    fi
    return $rc
}

overall=0
run_one fairnnv   FairNNV   matlab run_fairnnv.m         || overall=$?
# ProbVer uses a bash driver that runs each instance in its own MATLAB
# process — see ProbVer/run_probver.sh — so an OOM in one instance doesn't
# kill the whole verification run.
run_one probver   ProbVer   shell  run_probver.sh        || overall=$?
run_one gnnv      GNNV      matlab run_gnn_experiments.m || overall=$?
run_one videostar VideoStar matlab run_zoomin_4f.m       || overall=$?
run_one modelstar ModelStar matlab run_expt_for_compute.m || overall=$?
# ToolComparison defaults to smoke mode (~12 min) inside run_all.sh so the
# end-to-end suite stays under ~30 min. Override with:
#   TOOLCOMPARISON_MODE=full bash run_all.sh
# CPU-only; requires the AI Verification Library (AIVL) Support Package, which
# both Docker flows auto-install via `mpm`
# (`Deep_Learning_Toolbox_Verification_Library`). The build host's MATLAB
# licence must entitle the Verification Library.
export TOOLCOMPARISON_MODE="${TOOLCOMPARISON_MODE:-smoke}"
run_one toolcomparison ToolComparison matlab run_toolcomparison.m || overall=$?

# Consolidate experiment-specific result artifacts into repeatability_logs/results/
# so a single bind mount of repeatability_logs captures everything a reviewer
# needs to cross-check against the paper tables.
mkdir -p "${LOG_DIR}/results"
for exp_dir in FairNNV ProbVer GNNV VideoStar ModelStar; do
    src="${SCRIPT_DIR}/${exp_dir}/results"
    if [ ! -d "$src" ]; then continue; fi
    # Prefer the most recently modified timestamped subdir (e.g. results/<ts>/),
    # which is where FairNNV / GNNV / VideoStar / ProbVer write. ModelStar
    # writes directly into results/ without a subdir; the find branch catches
    # both shapes.
    latest=$(ls -dt "$src"/*/ 2>/dev/null | head -1)
    if [ -n "$latest" ]; then
        cp -r "$latest" "${LOG_DIR}/results/${exp_dir}" 2>/dev/null || true
    else
        cp -r "$src" "${LOG_DIR}/results/${exp_dir}" 2>/dev/null || true
    fi
done
# ToolComparison: not timestamped; .mat files in results/, rendered tables in
# tables/out/. Copy both so reviewers can compare table_main.{tex,txt} against
# the paper's aivl_comparison.tex.
if [ -d "${SCRIPT_DIR}/ToolComparison/results" ]; then
    mkdir -p "${LOG_DIR}/results/ToolComparison"
    cp -r "${SCRIPT_DIR}/ToolComparison/results/." \
        "${LOG_DIR}/results/ToolComparison/" 2>/dev/null || true
fi
if [ -d "${SCRIPT_DIR}/ToolComparison/tables/out" ]; then
    mkdir -p "${LOG_DIR}/results/ToolComparison/tables"
    cp -r "${SCRIPT_DIR}/ToolComparison/tables/out/." \
        "${LOG_DIR}/results/ToolComparison/tables/" 2>/dev/null || true
fi

{
    echo
    echo "=========================== SUMMARY ============================"
    column -t -s, "$SUMMARY_CSV" 2>/dev/null || cat "$SUMMARY_CSV"
    echo "================================================================"
    echo "Outputs written to repeatability_logs/ inside the container:"
    echo "    run.log                       consolidated, filtered terminal log"
    echo "    summary.csv                   per-experiment status + wall time"
    echo "    results/<experiment>/         per-experiment CSV / .mat / .tex artifacts"
    echo ""
    echo "If you followed the README docker-run command"
    echo "(-v \"\$PWD/results\":/out, with the trailing cp), these are now"
    echo "available on your HOST at:"
    echo "    \$PWD/results/repeatability_logs/"
} | tee -a "$RUN_LOG"
exit "$overall"
