#!/bin/bash
# run_all.sh — autonomous R2025b full ToolComparison sweep, NNV + AIVL.
#
# Runs every benchmark in the ToolComparison grid inside the nnv3.0:r2025b
# Docker image with the Vanderbilt MATLAB license. Sequential, non-interactive,
# resume-safe.
#
# Per benchmark:
#   1. matlab -batch inside container with bind-mounted host repo
#   2. NNV patches applied (idempotent)
#   3. AIVL preflight on the FIRST benchmark only (verifies
#      verifyNetworkRobustness actually runs, not just resolves)
#   4. run_toolcomparison('mode','default','tools',{'nnv','aivl'},'benchmarks',{<bench>})
#   5. Parse log for status counts and known issues
#   6. Append summary + issues to ISSUES.md
#
# Usage:
#   ./run_all.sh                   # all 6, nnv + aivl
#   ./run_all.sh nnv               # nnv only, all 6
#   ./run_all.sh nnv,aivl acas_xu_p3 rl
#
# Assumes the image `nnv3.0:r2025b` exists (run build_image.sh first).

set -uo pipefail

IMAGE_TAG="nnv3.0:r2025b"
LICENSE_SERVER="27009@licenseserver.it.vanderbilt.edu"
HOST_REPO="/home/verivital/Anne/dev/nnv3"
GUEST_REPO="/home/matlab/nnv"
TC_REL="code/nnv/examples/NNV3.0/ToolComparison"
TC_HOST="${HOST_REPO}/${TC_REL}"
TC_GUEST="${GUEST_REPO}/${TC_REL}"

LOG_DIR="$TC_HOST/logs"
ISSUES="$TC_HOST/ISSUES.md"
mkdir -p "$LOG_DIR"

# Image must exist.
if ! docker image inspect "$IMAGE_TAG" >/dev/null 2>&1; then
    echo "ERROR: image '$IMAGE_TAG' not found. Run build_image.sh first."
    exit 1
fi

# Args.
TOOLS="${1:-nnv,aivl}"
DEFAULT_BENCHES=(acas_xu_p3 acas_xu_p4 rl oval21 collins_rul mnist_resnet8)
if [[ $# -ge 2 ]]; then
    shift
    BENCHES=("$@")
else
    BENCHES=("${DEFAULT_BENCHES[@]}")
fi

TOOLS_CELL=$(echo "$TOOLS" | sed "s/,/','/g")
TOOLS_CELL="{'${TOOLS_CELL}'}"

# Per-benchmark wall budget (bash-side; per-instance budget is in instances.csv)
declare -A WALL
WALL[acas_xu_p3]=7200       # 2 h — exact-star dominates (300s/inst × 20)
WALL[acas_xu_p4]=7200       # 2 h — same profile as p3
WALL[rl]=3600               # 1 h — fast FC, sub-second per inst typical (50 inst now)
WALL[oval21]=7200           # 2 h — 30 inst × 4 algos; CIFAR ResNet, area variants
WALL[collins_rul]=7200      # 2 h — 62 inst × 4 algos; mostly fast under restricted grid
WALL[mnist_resnet8]=3600    # 1 h — 100 calls (25 imgs × 4 eps); inline, no parpool overhead

# Initialize ISSUES.md header if absent.
if [[ ! -f "$ISSUES" ]]; then
    cat > "$ISSUES" <<EOF
# ToolComparison R2025b default run notes

Run via: $0 $*
Image:   $IMAGE_TAG
License: $LICENSE_SERVER
Started: $(date -Iseconds)

EOF
fi

run_in_container() {
    # $1: matlab batch command string
    docker run --rm \
        -v "${HOST_REPO}:${GUEST_REPO}" \
        -e "MLM_LICENSE_FILE=${LICENSE_SERVER}" \
        -w "${TC_GUEST}" \
        --user matlab \
        "$IMAGE_TAG" \
        bash -lc "matlab -batch \"$1\""
}

# ---- AIVL preflight (runs once, before benchmark loop) ----
echo
echo "═══════════════════════════════════════════════════════════════════"
echo " AIVL preflight — verify verifyNetworkRobustness runs on R2025b"
echo "═══════════════════════════════════════════════════════════════════"
PREFLIGHT_LOG="$LOG_DIR/_preflight_aivl.log"
PREFLIGHT_CMD="cd ${TC_GUEST}/utils; addpath_shared; layers = [featureInputLayer(3); fullyConnectedLayer(2)]; net = dlnetwork(layers); xl = dlarray(single([-1;-1;-1]),'CB'); xu = dlarray(single([1;1;1]),'CB'); try; r = verifyNetworkRobustness(net, xl, xu, 1); fprintf('AIVL_OK %s', string(r)); catch ME; fprintf('AIVL_FAIL %s', ME.message); end"
docker run --rm \
    -v "${HOST_REPO}:${GUEST_REPO}" \
    -e "MLM_LICENSE_FILE=${LICENSE_SERVER}" \
    -w "${TC_GUEST}" \
    --user matlab \
    "$IMAGE_TAG" \
    bash -lc "matlab -batch \"${PREFLIGHT_CMD}\"" > "$PREFLIGHT_LOG" 2>&1

if grep -q "AIVL_OK" "$PREFLIGHT_LOG"; then
    echo "  ✓ AIVL preflight passed: $(grep AIVL_OK $PREFLIGHT_LOG | head -1)"
    AIVL_OK=1
else
    AIVL_OK=0
    echo "  ✗ AIVL preflight FAILED (see logs/_preflight_aivl.log)"
    grep -E "AIVL_FAIL|Error" "$PREFLIGHT_LOG" | head -5
    if [[ "$TOOLS" == *"aivl"* ]]; then
        echo
        echo "AIVL is in your tools list but doesn't run. Dropping AIVL from the tools."
        TOOLS=$(echo "$TOOLS" | sed 's/,aivl//; s/aivl,//; s/^aivl$/nnv/')
        TOOLS_CELL=$(echo "$TOOLS" | sed "s/,/','/g")
        TOOLS_CELL="{'${TOOLS_CELL}'}"
    fi
fi

{
    echo
    echo "## AIVL preflight"
    echo
    echo "- Result: $(grep -E 'AIVL_OK|AIVL_FAIL' $PREFLIGHT_LOG | head -1)"
    echo "- Tools after preflight: $TOOLS"
    echo "- Log: \`logs/_preflight_aivl.log\`"
} >> "$ISSUES"

# ---- Benchmark loop ----
TOTAL=${#BENCHES[@]}
IDX=0
GLOBAL_T0=$(date +%s)

for bench in "${BENCHES[@]}"; do
    IDX=$((IDX+1))
    wall_cap=${WALL[$bench]:-7200}
    cat <<EOF

═══════════════════════════════════════════════════════════════════
 Benchmark $IDX / $TOTAL : $bench
   tools:      $TOOLS
   wall cap:   ${wall_cap}s ($((wall_cap/60)) min)
   log file:   logs/${bench}.log
   results:    results/${bench}.mat
═══════════════════════════════════════════════════════════════════
EOF

    START=$(date +%s)
    CMD="addpath('${TC_GUEST}/utils'); addpath_shared; u = tool_utils; n = u.purge_status('${TC_GUEST}/results/${bench}.mat','error'); fprintf('purged %d error rows\\\\n', n); cd ${GUEST_REPO}/code/nnv; evalc('startup_nnv'); cd ${TC_GUEST}; run_toolcomparison('mode','default','tools',${TOOLS_CELL},'benchmarks',{'${bench}'})"

    timeout ${wall_cap}s docker run --rm \
        -v "${HOST_REPO}:${GUEST_REPO}" \
        -e "MLM_LICENSE_FILE=${LICENSE_SERVER}" \
        -w "${TC_GUEST}" \
        --user matlab \
        "$IMAGE_TAG" \
        bash -lc "matlab -batch \"${CMD}\"" > "$LOG_DIR/${bench}.log" 2>&1
    RC=$?
    END=$(date +%s)
    WALL_S=$((END - START))

    # Status counts via grep. Verdict appears at end-of-line for fast algos
    # (one-line log entry) and at line-start for exact-star (parpool spawn
    # message splits the line). Match the verdict anywhere on a line ending
    # with the timing parenthetical, e.g. "verified (0.48 s)" or
    # "timeout (>300 s)" — using word boundary to avoid false matches.
    nV=$(grep -cE '\<verified \([0-9.>]+ s\)' "$LOG_DIR/${bench}.log" 2>/dev/null || true)
    nX=$(grep -cE '\<violated \([0-9.>]+ s\)' "$LOG_DIR/${bench}.log" 2>/dev/null || true)
    nU=$(grep -cE '\<unknown \([0-9.>]+ s\)'  "$LOG_DIR/${bench}.log" 2>/dev/null || true)
    nT=$(grep -cE '\<timeout \([0-9.>]+ s\)'  "$LOG_DIR/${bench}.log" 2>/dev/null || true)
    nE=$(grep -cE '\<error \('                "$LOG_DIR/${bench}.log" 2>/dev/null || true)
    nTotal=$((nV + nX + nU + nT + nE))

    echo "  --- $bench done (rc=$RC, wall=${WALL_S}s = $((WALL_S/60))m$((WALL_S%60))s) ---"
    echo "  Instances: $nTotal | V=$nV X=$nX ?=$nU T/O=$nT Err=$nE"

    # Append to ISSUES.md
    {
        echo
        echo "## $bench"
        echo
        echo "- Wall: ${WALL_S}s ($(echo "scale=1; $WALL_S/60" | bc) min)"
        echo "- Container exit: $RC"
        echo "- Instances logged: $nTotal | V=$nV X=$nX ?=$nU T/O=$nT Err=$nE"
        echo "- Log: \`logs/${bench}.log\`"
        if grep -qE "Error using|Unrecognized|Unsupported|Unable to communicate|cannot be started from a worker|LP solver error" "$LOG_DIR/${bench}.log"; then
            echo "- ⚠ Issues detected:"
            grep -E "Error using|Unrecognized|Unsupported|Unable to communicate|cannot be started from a worker|LP solver error" \
                "$LOG_DIR/${bench}.log" | sort -u | head -5 | sed 's/^/  - /'
        fi
    } >> "$ISSUES"
done

GLOBAL_END=$(date +%s)
GLOBAL_WALL=$((GLOBAL_END - GLOBAL_T0))

cat <<EOF

═══════════════════════════════════════════════════════════════════
 All done. Total wall: ${GLOBAL_WALL}s ($((GLOBAL_WALL/60)) min)
   Per-bench logs: logs/
   Per-bench mats: results/
   Issues log:     ISSUES.md
═══════════════════════════════════════════════════════════════════

Rendering consolidated table…
EOF

run_in_container "addpath('${TC_GUEST}/utils'); addpath_shared; cd ${TC_GUEST}/tables; make_table_main" 2>&1 \
    | tail -25 > "$LOG_DIR/_make_table.log"
cat "$LOG_DIR/_make_table.log" | tail -25
