#!/usr/bin/env bash
# NNV3.0 — top-level orchestrator for the paper repeatability artifact.
#
# Drives the five example folders (FairNNV, ProbVer, VideoStar, GNNV,
# ModelStar) through a uniform interface. Each example continues to use
# its own internal style (see per-folder README), so this script is a
# thin shell over five matlab -batch invocations.
#
# Usage:
#   ./run_all.sh                   # smoke (default), all examples that don't need a GPU
#   ./run_all.sh --smoke           # alias for the above
#   ./run_all.sh --full            # paper-faithful reproduction (long; hours)
#   ./run_all.sh --no-gpu          # skip ProbVer (which requires a GPU)
#   ./run_all.sh --only gnnv       # run a single example (gnnv|fairnnv|probver|videostar|modelstar)
#   ./run_all.sh --smoke --no-gpu  # combinations are allowed
#
# Behavior:
#   - Each example is wrapped so a failure does NOT abort the rest.
#   - Per-example logs go to results/run_all_<timestamp>/<example>.log
#   - Pre-flight checks (license reachable, free disk, GPU) run once up front.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ----------------------- Argument parsing -----------------------

mode="smoke"
no_gpu=0
only=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --smoke) mode="smoke"; shift ;;
        --full)  mode="full"; shift ;;
        --no-gpu) no_gpu=1; shift ;;
        --only)
            shift
            only="${1:-}"
            if [[ -z "$only" ]]; then
                echo "ERROR: --only requires a value (gnnv|fairnnv|probver|videostar|modelstar)" >&2
                exit 2
            fi
            shift
            ;;
        -h|--help)
            sed -n '2,21p' "$0"
            exit 0
            ;;
        *) echo "ERROR: unknown flag $1" >&2; exit 2 ;;
    esac
done

valid_examples=(gnnv fairnnv probver videostar modelstar)
if [[ -n "$only" ]]; then
    if ! printf '%s\n' "${valid_examples[@]}" | grep -qx "$only"; then
        echo "ERROR: --only must be one of: ${valid_examples[*]}" >&2
        exit 2
    fi
fi

# ----------------------- Pre-flight checks -----------------------

ts="$(date +%y%m%d-%H%M%S)"
out_root="$SCRIPT_DIR/results/run_all_${ts}"
mkdir -p "$out_root"

echo "================================================================"
echo "  NNV3.0 — run_all  (mode=${mode}, no_gpu=${no_gpu}, only=${only:-all})"
echo "================================================================"
echo "Output root: $out_root"
echo

# License: presence of MLM_LICENSE_FILE OR a license server reachable on tcp/27009.
license_warn=""
if [[ -z "${MLM_LICENSE_FILE:-}" ]]; then
    license_warn="MLM_LICENSE_FILE is not set; relying on default MATLAB license setup."
fi
[[ -n "$license_warn" ]] && echo "WARN: $license_warn"

# Free disk: need at least 5 GB under SCRIPT_DIR (Docker results + decompressed assets).
free_kb=$(df -k "$SCRIPT_DIR" | awk 'NR==2{print $4}')
if [[ "${free_kb:-0}" -lt 5000000 ]]; then
    echo "WARN: less than 5 GB free under $SCRIPT_DIR (current: $((free_kb/1024)) MB)."
fi

# GPU presence: nvidia-smi must exist and report at least one GPU.
has_gpu=0
if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L 2>/dev/null | grep -q GPU; then
    has_gpu=1
fi
if [[ "$no_gpu" -eq 0 && "$has_gpu" -eq 0 ]]; then
    echo "WARN: no GPU detected (nvidia-smi). ProbVer will fail without one. Pass --no-gpu to skip."
fi

echo

# ----------------------- Per-example invocations -----------------------

# matlab_run <log-file> <matlab-command>
# Returns 0 on success, non-zero on failure. Never aborts the script.
matlab_run() {
    local log="$1"; shift
    local cmd="$1"; shift
    echo ">> matlab -batch \"$cmd\""
    if matlab -nodisplay -batch "$cmd" >"$log" 2>&1; then
        echo "   OK"
        return 0
    else
        echo "   FAILED — see $log"
        return 1
    fi
}

want() {
    # Should we run this example?
    local name="$1"
    if [[ -n "$only" && "$only" != "$name" ]]; then return 1; fi
    return 0
}

# Smoke = quick sanity subset (~minutes). Full = paper-faithful sweep.
# Each per-folder runner has a config-guard, so overrides set BEFORE
# `run('...')` are honored. Full uses each script's default config.
fairnnv_smoke="config.modelList = {'AC-1'}; config.numObs = 10; config.epsilon_individual = [0.01]; config.timeout = 120;"
fairnnv_full=""    # use script's defaults
probver_smoke="numSamples = 1; nRand = 50;"
probver_full=""    # numSamples = 3 by default
videostar_smoke="config.sampleIndices = 1:1; config.epsilon = [1/255]; config.timeout = 300;"
videostar_full=""  # sampleIndices = 1:10, three epsilons
gnnv_smoke="run_gnn_experiments('num_graphs', 3, 'architectures', {'gcn'}, 'node_epsilons', [1e-3]);"
gnnv_full="run_gnn_experiments();"
modelstar_smoke="n_layers_to_run_for_from_yaml_file = 1;"
modelstar_full=""    # use script's defaults (3 layers: fc_4/fc_5/fc_6)

run_fairnnv() {
    want fairnnv || return 0
    local log="$out_root/fairnnv.log"
    echo "[FairNNV]"
    local cfg="$fairnnv_smoke"
    [[ "$mode" == "full" ]] && cfg="$fairnnv_full"
    # Set overrides BEFORE the runner so its config-guard picks them up.
    matlab_run "$log" "cd('$SCRIPT_DIR/FairNNV'); $cfg run('run_fairnnv.m');"
}

run_probver() {
    want probver || return 0
    if [[ "$no_gpu" -eq 1 ]]; then
        echo "[ProbVer] SKIPPED: --no-gpu set (cp-star is much slower on CPU)"
        return 0
    fi
    local log="$out_root/probver.log"
    echo "[ProbVer]"
    local cfg="$probver_smoke"
    [[ "$mode" == "full" ]] && cfg="$probver_full"
    matlab_run "$log" "cd('$SCRIPT_DIR/ProbVer'); $cfg run('run_probver.m');"
}

run_videostar() {
    want videostar || return 0
    local log="$out_root/videostar.log"
    echo "[VideoStar]"
    local cfg="$videostar_smoke"
    [[ "$mode" == "full" ]] && cfg="$videostar_full"
    matlab_run "$log" "cd('$SCRIPT_DIR/VideoStar'); $cfg run('run_zoomin_4f.m');"
}

run_gnnv() {
    want gnnv || return 0
    local log="$out_root/gnnv.log"
    echo "[GNNV]"
    local cmd="$gnnv_smoke"
    [[ "$mode" == "full" ]] && cmd="$gnnv_full"
    matlab_run "$log" "cd('$SCRIPT_DIR/GNNV'); $cmd"
}

run_modelstar() {
    want modelstar || return 0
    local log="$out_root/modelstar.log"
    echo "[ModelStar]"
    local cfg="$modelstar_smoke"
    [[ "$mode" == "full" ]] && cfg="$modelstar_full"
    matlab_run "$log" "cd('$SCRIPT_DIR/ModelStar'); $cfg run('run_expt_for_compute.m');"
}

# Order chosen so the cheapest examples (GNNV, FairNNV, ModelStar) run first,
# giving fast feedback before the slower ones.
overall_rc=0
run_gnnv      || overall_rc=1
run_fairnnv   || overall_rc=1
run_modelstar || overall_rc=1
run_videostar || overall_rc=1
run_probver   || overall_rc=1

echo
echo "================================================================"
if [[ "$overall_rc" -eq 0 ]]; then
    echo "  All requested examples completed."
else
    echo "  At least one example FAILED. See per-example logs in:"
    echo "  $out_root"
fi
echo "================================================================"
exit "$overall_rc"
