#!/bin/bash
# run_step_by_step.sh — execute the default-mode grid one benchmark at a
# time, pausing between each so the user can inspect results and log issues.
#
# Usage:
#   ./run_step_by_step.sh                  # all 6 benchmarks, nnv+aivl
#   ./run_step_by_step.sh nnv              # nnv only
#   ./run_step_by_step.sh nnv acas_xu_p3 rl  # nnv only, specific
#
# Per benchmark:
#   1. Print scope + ETA
#   2. matlab -batch run_toolcomparison('benchmarks',{<bench>},'mode','default',…)
#   3. Logs to logs/<bench>.log
#   4. Parse log for status counts (no extra matlab launch)
#   5. Append summary + auto-detected issues to ISSUES.md
#   6. Pause for user (Enter=next, q=quit)
#
# Only one MATLAB license checkout per benchmark.

set -uo pipefail

TC_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$TC_DIR" && git rev-parse --show-toplevel)"
LOG_DIR="$TC_DIR/logs"
ISSUES="$TC_DIR/ISSUES.md"
mkdir -p "$LOG_DIR"

# -------- arg parsing --------
TOOLS="nnv,aivl"
DEFAULT_BENCHES=(acas_xu_p3 acas_xu_p4 rl oval21 collins_rul mnist_resnet8)
ARGS=("$@")
BENCHES=()

if [[ ${#ARGS[@]} -ge 1 ]] && [[ "${ARGS[0]}" =~ ^(nnv|aivl|nnv,aivl|aivl,nnv)$ ]]; then
    TOOLS="${ARGS[0]}"
    ARGS=("${ARGS[@]:1}")
fi
if [[ ${#ARGS[@]} -gt 0 ]]; then
    BENCHES=("${ARGS[@]}")
else
    BENCHES=("${DEFAULT_BENCHES[@]}")
fi

TOOLS_CELL=$(echo "$TOOLS" | sed "s/,/','/g")
TOOLS_CELL="{'${TOOLS_CELL}'}"

# Per-benchmark ETA hints
declare -A ETA
ETA[acas_xu_p3]="20-25 min  (20 inst × 3 NNV + 1 AIVL; exact-star at ~30 s/inst dominates)"
ETA[acas_xu_p4]="20-25 min  (20 inst × 3 NNV + 1 AIVL; same profile as p3)"
ETA[rl]="10-15 min  (50 inst × 3 NNV + 1 AIVL; fast FC, sub-second per inst)"
ETA[oval21]="40-60 min  (30 inst × 2 NNV + 1 AIVL; CIFAR ResNet, area variants)"
ETA[collins_rul]="20-40 min  (62 inst × 2 NNV + 1 AIVL; mostly fast)"
ETA[mnist_resnet8]="30-40 min  (100 calls × ~8 s NNV + ~1.5 s AIVL per call)"

if [[ ! -f "$ISSUES" ]]; then
    cat > "$ISSUES" <<EOF
# ToolComparison default-run notes

Tracks per-benchmark wall time, status mix, and any issues observed
during the step-by-step default run.

EOF
fi

TOTAL=${#BENCHES[@]}
IDX=0
GLOBAL_T0=$(date +%s)

for bench in "${BENCHES[@]}"; do
    IDX=$((IDX+1))
    eta="${ETA[$bench]:-unknown}"
    cat <<EOF

═══════════════════════════════════════════════════════════════════════════
 Benchmark $IDX / $TOTAL : $bench
   tools:      $TOOLS
   ETA:        $eta
   log file:   logs/${bench}.log
   results:    results/${bench}.mat
═══════════════════════════════════════════════════════════════════════════
EOF

    START=$(date +%s)
    cd "$REPO_ROOT"

    matlab -batch "cd code/nnv; evalc('startup_nnv'); cd $TC_DIR; run_toolcomparison('mode','default','tools',${TOOLS_CELL},'benchmarks',{'$bench'})" \
        > "$LOG_DIR/${bench}.log" 2>&1
    RC=$?
    END=$(date +%s)
    WALL=$((END - START))

    # Status counts via grep. Matlab interleaves the "[bench tool alg]…" prefix
    # with parpool-spawn messages so the verdict word ends up on its own line:
    #   "verified (7.7 s)"
    #   "timeout (>30 s)"
    nV=$(grep -cE '^verified \('  "$LOG_DIR/${bench}.log" || true)
    nX=$(grep -cE '^violated \('  "$LOG_DIR/${bench}.log" || true)
    nU=$(grep -cE '^unknown \('   "$LOG_DIR/${bench}.log" || true)
    nT=$(grep -cE '^timeout \('   "$LOG_DIR/${bench}.log" || true)
    nE=$(grep -cE '^error \('     "$LOG_DIR/${bench}.log" || true)
    nTotal=$((nV + nX + nU + nT + nE))
    SUMMARY="  Instances: ${nTotal} total | verified=${nV} violated=${nX} unknown=${nU} timeout=${nT} error=${nE}"

    echo
    echo "  --- $bench done (exit=$RC, wall=${WALL}s = $((WALL/60))m$((WALL%60))s) ---"
    echo "$SUMMARY"

    # Append to ISSUES.md
    {
        echo
        echo "## $bench"
        echo
        echo "- Wall: ${WALL}s ($(echo "scale=1; $WALL/60" | bc) min)"
        echo "- Exit code: $RC"
        echo "- $SUMMARY"
        echo "- Log: \`logs/${bench}.log\`"
        # Auto-detect common issues
        if grep -qE "Error using|Unrecognized|Unsupported|Unable to communicate" "$LOG_DIR/${bench}.log"; then
            echo "- ⚠ Issues detected:"
            grep -E "Error using|Unrecognized|Unsupported|Unable to communicate" "$LOG_DIR/${bench}.log" \
                | sort -u | head -3 | sed 's/^/  - /'
        fi
    } >> "$ISSUES"

    if [[ $IDX -lt $TOTAL ]]; then
        echo
        read -r -p "Press Enter to continue to next benchmark, [q] to quit: " key
        case "$key" in
            q|Q) echo "User quit."; break ;;
        esac
    fi
done

GLOBAL_END=$(date +%s)
GLOBAL_WALL=$((GLOBAL_END - GLOBAL_T0))

cat <<EOF

═══════════════════════════════════════════════════════════════════════════
 All done. Total wall: ${GLOBAL_WALL}s ($((GLOBAL_WALL/60)) min)
   Per-benchmark logs: logs/
   Per-benchmark mats: results/
   Issues log:         ISSUES.md
═══════════════════════════════════════════════════════════════════════════

Rendering consolidated table…
EOF

cd "$REPO_ROOT"
matlab -batch "cd code/nnv; evalc('startup_nnv'); cd $TC_DIR/tables; addpath('$TC_DIR/utils'); make_table_main" 2>&1 | tail -25
