#!/usr/bin/env bash
# run_probver.sh - Memory-isolated driver for the yolo_2023 ProbVer suite.
#
# Each TinyYOLO instance runs in its own MATLAB process (one MATLAB ->
# verify_one_instance(idx, ...) -> exit). When `Prob_reach` blows past
# host memory and the kernel SIGKILLs MATLAB, only that instance is
# lost; the next one starts fresh and the CSV keeps growing.
#
# Output: results_summary.csv next to this script. Header is written
# fresh; one row per instance, with status `oom` recorded for crashed
# instances so the suite always finishes with a complete row count.
#
# Tunables (env vars):
#   PROBVER_NUM_SAMPLES  number of instances to verify (default 3)
#   PROBVER_SEED         seed for instance selection      (default 42)
#   PROBVER_NRAND        falsification samples per inst.  (default 100)
#   MATLAB               matlab binary on PATH            (default 'matlab')

set -uo pipefail

NUM_SAMPLES="${PROBVER_NUM_SAMPLES:-3}"
SEED="${PROBVER_SEED:-42}"
NRAND="${PROBVER_NRAND:-100}"
MATLAB_BIN="${MATLAB:-matlab}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="${SCRIPT_DIR}/yolo_2023"
OUT_CSV="${SCRIPT_DIR}/results_summary.csv"

if [[ ! -d "$BENCH_DIR" ]]; then
    echo "[probver] benchmark dir not found: $BENCH_DIR" >&2
    exit 2
fi

echo "[probver] decompressing benchmark files (idempotent)..."
if [[ -f "$BENCH_DIR/onnx/TinyYOLO.onnx.gz" && ! -f "$BENCH_DIR/onnx/TinyYOLO.onnx" ]]; then
    gunzip -k "$BENCH_DIR/onnx/TinyYOLO.onnx.gz"
fi
shopt -s nullglob
for gz in "$BENCH_DIR/vnnlib/"*.vnnlib.gz; do
    plain="${gz%.gz}"
    if [[ ! -f "$plain" ]]; then
        gunzip -k "$gz"
    fi
done
shopt -u nullglob

# Pick instance indices in a single, cheap MATLAB session that allocates
# almost nothing. Output one line of the form `PICKED: 27 52 68`.
echo "[probver] picking ${NUM_SAMPLES} indices (seed=${SEED})..."
PICK_OUT="$( "$MATLAB_BIN" -nodisplay -batch "\
    fid = fopen('${BENCH_DIR}/instances.csv','r'); \
    n = 0; while ~feof(fid), line = fgetl(fid); if ischar(line) && ~isempty(strtrim(line)), n = n + 1; end; end; fclose(fid); \
    rng(${SEED}); \
    idx = sort(randperm(n, min(${NUM_SAMPLES}, n))); \
    fprintf('PICKED:'); fprintf(' %d', idx); fprintf('\\n');" 2>&1 \
    | tr -d '\r' | grep '^PICKED:' || true )"
INDICES="$(echo "$PICK_OUT" | sed 's/^PICKED://' | xargs)"

if [[ -z "$INDICES" ]]; then
    echo "[probver] could not pick instances; MATLAB output was:" >&2
    echo "$PICK_OUT" >&2
    exit 3
fi
echo "[probver] selected: $INDICES"

# Fresh CSV with header.
echo "index,onnx,vnnlib,status,time,error" > "$OUT_CSV"

overall=0
total_t0=$(date -u +%s)
for idx in $INDICES; do
    t0=$(date -u +%s)
    echo
    echo "[probver] === instance $idx ==="
    set +e
    "$MATLAB_BIN" -nodisplay -batch \
        "verify_one_instance(${idx}, '${BENCH_DIR}', '${OUT_CSV}', ${NRAND})"
    rc=$?
    set -e
    t1=$(date -u +%s)
    elapsed=$(( t1 - t0 ))

    if [[ $rc -ne 0 ]]; then
        # MATLAB crashed (often SIGKILL = 137 from the OOM killer). Append
        # a marker row so the CSV reflects exactly numSamples lines.
        instance_meta="$(awk -F, -v i="$idx" 'NR==i { gsub(/^[ \t]+|[ \t]+$/, "", $1); gsub(/^[ \t]+|[ \t]+$/, "", $2); print $1","$2 }' "$BENCH_DIR/instances.csv")"
        if [[ -z "$instance_meta" ]]; then
            instance_meta=",,"
        fi
        if [[ "$rc" -eq 137 ]]; then
            why="MATLAB SIGKILLed (likely OOM during reachability)"
            stat="oom"
        else
            why="MATLAB exit ${rc}"
            stat="error"
        fi
        printf '%s,%s,%s,%d,%s\n' "$idx" "$instance_meta" "$stat" "$elapsed" "$why" >> "$OUT_CSV"
        echo "[probver]   recorded as '$stat' (rc=$rc) after ${elapsed}s"
        # Do not abort — keep going through the remaining instances.
        overall=$(( overall + 1 ))
    else
        echo "[probver]   ok in ${elapsed}s"
    fi
done

total_elapsed=$(( $(date -u +%s) - total_t0 ))

echo
echo "===================== ProbVer SUMMARY ====================="
column -t -s, "$OUT_CSV" 2>/dev/null || cat "$OUT_CSV"
echo "==========================================================="
echo "Total wall-clock: ${total_elapsed}s   Failed instances: ${overall}"
exit "$overall"
