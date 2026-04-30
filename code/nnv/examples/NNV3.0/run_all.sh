#!/usr/bin/env bash
# Run the full NNV3.0 repeatability suite (FairNNV + ProbVer + GNNV +
# VideoStar + ModelStar) in one shot. Each experiment runs in its own
# MATLAB session so a crash in one does not lose the others.
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
#                                       'full' takes ~6 h and won't fit a
#                                       typical CodeOcean compute window.
#
# Exit code is 0 only if every selected experiment finished cleanly. The final
# summary table (and a CSV at $NNV3_LOG_DIR/summary.csv) lists wall-clock time
# per experiment.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NNV_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"  # .../code/nnv
LOG_DIR="${NNV3_LOG_DIR:-${SCRIPT_DIR}/repeatability_logs}"
MATLAB="${MATLAB:-matlab}"
SKIP="${NNV3_SKIP:-}"

mkdir -p "$LOG_DIR"
SUMMARY_CSV="${LOG_DIR}/summary.csv"
echo "experiment,status,wall_seconds,log" > "$SUMMARY_CSV"

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

# Forward-compat call is repeated in each script too; we set it here as well
# so MATLAB doesn't error on Blackwell / RTX 5090 hosts. The call is wrapped
# in try/catch and is a no-op on hosts without an NVIDIA GPU.
MATLAB_PRELUDE="addpath(genpath('${NNV_ROOT}')); try, parallel.gpu.enableCUDAForwardCompatibility(true); catch; end;"

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
        printf '\n=== %-12s SKIPPED (%s) ===\n' "$name" "$reason"
        echo "${name},skipped,0,${reason}" >> "$SUMMARY_CSV"
        return 0
    fi

    local logfile="${LOG_DIR}/${name}.log"
    printf '\n=== %-12s start: %s ===\n' "$name" "$(date -u +%FT%TZ)"
    local t0
    t0=$(date -u +%s)
    if [[ "$kind" == "matlab" ]]; then
        (
            cd "${SCRIPT_DIR}/${subdir}" || exit 99
            "$MATLAB" -nodisplay -batch "${MATLAB_PRELUDE} run('${entry}'); exit()"
        ) 2>&1 | tee "$logfile"
    else
        (
            cd "${SCRIPT_DIR}/${subdir}" || exit 99
            bash "${SCRIPT_DIR}/${subdir}/${entry}"
        ) 2>&1 | tee "$logfile"
    fi
    local rc=${PIPESTATUS[0]}
    local t1
    t1=$(date -u +%s)
    local elapsed=$(( t1 - t0 ))
    if [[ $rc -eq 0 ]]; then
        printf '=== %-12s OK in %ds (log: %s) ===\n' "$name" "$elapsed" "$logfile"
        echo "${name},ok,${elapsed},${logfile}" >> "$SUMMARY_CSV"
    else
        printf '=== %-12s FAILED (rc=%d) in %ds (log: %s) ===\n' "$name" "$rc" "$elapsed" "$logfile"
        echo "${name},failed,${elapsed},${logfile}" >> "$SUMMARY_CSV"
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
# CPU-only; requires the AI Verification Toolbox Support Package staged
# at /home/matlab/addons/atva26-aivl.tar.gz (or the other paths checked
# by ToolComparison/scripts/toolbox_install.m).
export TOOLCOMPARISON_MODE="${TOOLCOMPARISON_MODE:-smoke}"
run_one toolcomparison ToolComparison matlab run_toolcomparison.m || overall=$?

echo
echo "=========================== SUMMARY ============================"
column -t -s, "$SUMMARY_CSV" 2>/dev/null || cat "$SUMMARY_CSV"
echo "================================================================"
echo "Per-experiment logs: $LOG_DIR"
exit "$overall"
