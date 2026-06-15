#!/bin/bash
# Pre-sweep benchmark SMOKE TEST. Runs ONE representative instance per benchmark
# (run_all_benchmarks 'first' mode) with a short timeout, then reports which
# benchmarks ERROR (load / ONNX-import / dispatch failures) -- as distinct from
# legitimate timeouts/unknowns.
#
# WHY: a full sweep is hours of compute; a benchmark that errors on instance 1
# errors on all of them. Run this BEFORE every sweep, and after ANY benchmark
# change or NNV code change -- ESPECIALLY once GPU-BaB is routed into
# run_vnncomp_instance -- so breakage is caught in minutes, not after a wasted sweep.
# It exercises the exact same run_vnncomp_instance path the sweep uses, so whatever
# routing the sweep does (incl. future GPU-BaB) is covered automatically.
#
# Exit code: 0 = all benchmarks ran clean; 1 = one or more ERRORED (list printed);
#            2 = harness failure (no results produced).
# Usage: bash smoke_test_benchmarks.sh [bench_root] [timeout_s]
set -uo pipefail

ROOT="${1:-$HOME/vnncomp2026_benchmarks/benchmarks}"
TIMEOUT="${2:-30}"
SUBDIR="$(cd "$(dirname "$0")/.." && pwd)"      # .../Submission/VNN_COMP2026
NNVROOT="$(cd "$SUBDIR/../../.." && pwd)"        # .../code/nnv
MATLAB_BIN="$(command -v matlab || echo /usr/local/matlab/bin/matlab)"

echo "[smoke] root=$ROOT  timeout=${TIMEOUT}s  runner=$SUBDIR/run_all_benchmarks.m"
[ -d "$ROOT" ] || { echo "[smoke] FATAL: benchmark root not found: $ROOT"; exit 2; }

# Run one instance per benchmark through the maintained runner.
"$MATLAB_BIN" -batch "cd('$NNVROOT'); startup_nnv; addpath('$SUBDIR'); cd('$SUBDIR'); run_all_benchmarks('$ROOT', $TIMEOUT, {}, 'first');" \
  || { echo "[smoke] FATAL: matlab run_all_benchmarks failed"; exit 2; }

# Find the newest results CSV the runner wrote and report errored benchmarks.
CSV="$(ls -t "$SUBDIR"/results_*.csv 2>/dev/null | head -1)"
[ -n "$CSV" ] || { echo "[smoke] FATAL: no results_*.csv produced"; exit 2; }
echo "[smoke] results: $CSV"

python3 - "$CSV" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
err = [r for r in rows if r.get('status_str','').strip().lower() == 'error']
print(f"[smoke] {len(rows)} benchmarks tested, {len(err)} ERRORED")
for r in err:
    print(f"  ERROR  {r.get('subfolder','?'):<40} {r.get('error_message','')[:110]}")
if not err:
    print("[smoke] OK -- all benchmarks ran clean; safe to launch the full sweep.")
sys.exit(1 if err else 0)
PY
