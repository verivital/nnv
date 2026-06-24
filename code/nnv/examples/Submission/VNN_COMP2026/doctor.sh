#!/bin/bash
# VNN-COMP 2026 NNV preflight "doctor" -- verify the eval/sweep environment is COMPLETE before a run,
# so a missing dep / wrong python / unset asset fails LOUDLY here instead of silently degrading a whole
# sweep to `unknown` or, worse, FAIL-OPENING the authoritative witness gate (-> a spurious sat stands =
# a -150). READ-ONLY: it imports + version-checks, never installs. Exit 0 = ready, 1 = a FATAL gap.
#
#   bash doctor.sh [bench_root]        # bench_root default: ../../../../../../vnncomp2026_benchmarks/benchmarks
#
# What it checks (and who needs it):
#   - python with onnx+onnxruntime+vnnlib   -> the witness gate (BOTH the official path AND the sweep) +
#                                              the python falsifiers; vnnlib-less => gate fail-opens.
#   - python with the full onnx2nnv stack   -> prepare_instance.sh manifest generation (lsnc/traffic/...).
#   - matlab on PATH                         -> every run.
#   - matlab.engine                          -> the OFFICIAL execute.py path only (sweeps use matlab -batch).
#   - key harness files + benchmark assets present.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REQ="$(realpath -m "$HERE/../../../tools/onnx2nnv_python/requirements.txt" 2>/dev/null || echo "$HERE/../../../tools/onnx2nnv_python/requirements.txt")"
BENCH="${1:-$HERE/../../../../../../vnncomp2026_benchmarks/benchmarks}"
BENCH="$(realpath -m "$BENCH" 2>/dev/null || echo "$BENCH")"

FATAL=0; WARN=0
pass(){ printf '  \033[32mPASS\033[0m %s\n' "$1"; }
warn(){ printf '  \033[33mWARN\033[0m %s\n' "$1"; WARN=$((WARN+1)); }
fail(){ printf '  \033[31mFAIL\033[0m %s\n' "$1"; FATAL=$((FATAL+1)); }
hdr(){  printf '\n== %s ==\n' "$1"; }

# Resolve the FIRST python that imports a given module set (mirrors python_exe / i_ort_vnnlib_python:
# NNV_ORT_PYTHON -> ~/taylor_venv -> python3 -> python). Echoes the interpreter or '' .
resolve_py(){
  local mods="$1" c
  for c in "${NNV_ORT_PYTHON:-}" "$HOME/taylor_venv/bin/python" python3 python; do
    [ -z "$c" ] && continue
    if "$c" -c "import $mods" >/dev/null 2>&1; then echo "$c"; return 0; fi
  done
  echo ""
}

hdr "VNN-COMP 2026 NNV environment doctor"
echo "  submission dir: $HERE"
echo "  bench root:     $BENCH"
echo "  requirements:   $REQ"

hdr "1) Witness-gate / falsifier python (onnx + onnxruntime + vnnlib)"
GPY="$(resolve_py 'onnx, onnxruntime, vnnlib')"
if [ -n "$GPY" ]; then
  pass "found: $GPY"
  V=$("$GPY" -c 'import onnx,onnxruntime as o; print(onnx.__version__, o.__version__)' 2>/dev/null)
  echo "        onnx=$(echo "$V" | awk '{print $1}')  onnxruntime=$(echo "$V" | awk '{print $2}')"
  # Compare to the requirements.txt pins; a divergence is non-fatal (the gate uses a lenient tol) but
  # is the class that bit adaptive_cruise (dev ort != witness-generation ort) -> surface it.
  if [ -f "$REQ" ]; then
    PIN_ONNX=$(grep -oE 'onnx==[0-9.]+' "$REQ" | cut -d= -f3)
    PIN_ORT=$(grep -oE 'onnxruntime==[0-9.]+' "$REQ" | cut -d= -f3)
    HAVE_ONNX=$(echo "$V" | awk '{print $1}'); HAVE_ORT=$(echo "$V" | awk '{print $2}')
    [ -n "$PIN_ONNX" ] && [ "$PIN_ONNX" != "$HAVE_ONNX" ] && warn "onnx $HAVE_ONNX != requirements pin $PIN_ONNX (eval box should match; dev tolerable)"
    [ -n "$PIN_ORT" ]  && [ "$PIN_ORT"  != "$HAVE_ORT"  ] && warn "onnxruntime $HAVE_ORT != requirements pin $PIN_ORT (the adaptive_cruise mismatch class)"
  fi
else
  fail "no python imports onnx+onnxruntime+vnnlib -> the witness gate FAIL-OPENS (spurious sats stand). Install: <python> -m pip install -r $REQ"
fi

hdr "2) Manifest importer python (onnx2nnv full stack)"
IPY="$(resolve_py 'numpy, scipy, onnx, onnxruntime, onnxsim, onnxoptimizer, vnnlib')"
if [ -n "$IPY" ]; then
  pass "found: $IPY"
else
  fail "no python imports the full onnx2nnv stack (numpy scipy onnx onnxruntime onnxsim onnxoptimizer vnnlib) -> manifest categories (lsnc/traffic/cgan/soundnessbench/nn4sys/vit) error. Install: <python> -m pip install -r $REQ"
fi

hdr "3) MATLAB"
if command -v matlab >/dev/null 2>&1; then
  pass "matlab on PATH: $(command -v matlab)"
else
  fail "matlab not on PATH (set it or pass MATLAB explicitly to the sweep launcher)"
fi

hdr "4) matlab.engine (official execute.py path only; sweeps use matlab -batch)"
EPY="$(resolve_py 'matlab.engine')"
if [ -n "$EPY" ]; then pass "matlab.engine importable via: $EPY"
else warn "matlab.engine not importable (FATAL only for the official execute.py path; dev sweeps via sweep_lambda are unaffected)"; fi

hdr "5) Harness files present"
for f in run_vnncomp_instance.m run_all_benchmarks.m authoritative_witness_gate.m \
         validate_witness_authoritative.py execute.py; do
  [ -f "$HERE/$f" ] && pass "$f" || fail "missing: $f"
done

hdr "6) Benchmark assets"
if [ -d "$BENCH" ]; then
  pass "bench root exists"
  N=$(find "$BENCH" -maxdepth 3 -name instances.csv 2>/dev/null | wc -l)
  [ "$N" -gt 0 ] && pass "instances.csv files found: $N" || warn "no instances.csv under bench root (run setup.sh + gunzip?)"
else
  warn "bench root not found: $BENCH (pass it as arg 1, or run setup.sh)"
fi

hdr "Summary"
printf '  %d fatal, %d warn\n' "$FATAL" "$WARN"
if [ "$FATAL" -gt 0 ]; then
  printf '\033[31mNOT READY\033[0m -- resolve the FAIL items above before running.\n'; exit 1
fi
printf '\033[32mREADY\033[0m -- environment looks complete%s.\n' "$([ "$WARN" -gt 0 ] && echo " (review WARNs)")"
exit 0
