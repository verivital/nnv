#!/bin/bash
# Automated fresh-AWS-box bootstrap for NNV VNN-COMP 2026 work. IDEMPOTENT (re-run safe).
#
# Pulls EVERYTHING from public repositories -- no box-to-box transfer. The ONLY thing
# provisioned out of band is the MATLAB license (it is a SECRET; never put it in a repo).
#
# Target: Ubuntu + MATLAB R2026a (the MathWorks AWS AMI, or a self-installed MATLAB on PATH).
# Run as the MATLAB user (e.g. ubuntu). Disk: budget ~60 GB -- the benchmark repos clone
# at ~2 GB each and decompress to ~13 GB (2025) + ~25 GB (2026).
#
# Usage:   bash setup_aws_box.sh
# Tunables (env): NNV_BRANCH (default the 2026 migration branch -- switch to master after merge),
#                 RUN_BENCH_SETUP=0 to skip the (slow) benchmark decompress, *_REPO overrides.
set -uo pipefail

# ---- config (override via env) ------------------------------------------------------------
NNV_REPO="${NNV_REPO:-https://github.com/ttj/nnv.git}"
NNV_BRANCH="${NNV_BRANCH:-master}"   # has VNN_COMP2026 self-contained submission
BENCH2025_REPO="${BENCH2025_REPO:-https://github.com/VNN-COMP/vnncomp2025_benchmarks.git}"
BENCH2026_REPO="${BENCH2026_REPO:-https://github.com/VNN-COMP/vnncomp2026_benchmarks.git}"
RESULTS2025_REPO="${RESULTS2025_REPO:-https://github.com/VNN-COMP/vnncomp2025_results.git}"
STATUS_REPO="${STATUS_REPO:-https://github.com/ttj/nnv-vnncomp2026-status.git}"
RUN_BENCH_SETUP="${RUN_BENCH_SETUP:-1}"                  # run each benchmark repo's own setup.sh

log(){ echo "[setup $(date -u +%H:%M:%S)] $*"; }
clone_or_pull(){  # $1 url  $2 dir  [$3 branch]
  if [ -d "$2/.git" ]; then
    log "refresh $(basename "$2")"; git -C "$2" fetch -q origin || true
    if [ -n "${3:-}" ]; then
      git -C "$2" checkout -q -f "$3" 2>/dev/null && git -C "$2" pull -q --ff-only origin "$3" 2>/dev/null || true
    fi
  else
    log "clone $(basename "$2") <- $1 ${3:+($3)}"; git clone -q ${3:+-b "$3"} "$1" "$2"
  fi
}

# ---- 0. system deps -----------------------------------------------------------------------
log "apt deps (git, python3-pip, unzip, wget)"
sudo apt-get update -y -q && sudo apt-get install -y -q git python3 python3-pip unzip wget

# ---- 1. NNV (the migration branch carries the self-contained VNN_COMP2026 submission) ------
clone_or_pull "$NNV_REPO" "$HOME/nnv" "$NNV_BRANCH"

# ---- 2. public data: benchmarks (2025 + 2026), official 2025 results, status repo ----------
clone_or_pull "$BENCH2025_REPO"   "$HOME/vnncomp2025_benchmarks"
clone_or_pull "$BENCH2026_REPO"   "$HOME/vnncomp2026_benchmarks"
clone_or_pull "$RESULTS2025_REPO" "$HOME/vnncomp2025_results"
clone_or_pull "$STATUS_REPO"      "$HOME/nnv-vnncomp2026-status"

# ---- 2b. each benchmark repo ships a setup.sh that decompresses/prepares its models --------
if [ "$RUN_BENCH_SETUP" = "1" ]; then
  for d in vnncomp2025_benchmarks vnncomp2026_benchmarks; do
    if [ -f "$HOME/$d/setup.sh" ]; then
      log "running $d/setup.sh (decompress -- slow)"; ( cd "$HOME/$d" && bash setup.sh ) || log "WARN: $d setup.sh returned non-zero"
    fi
  done
fi

# ---- 3. MATLAB license (SECRET; out of band, NEVER committed) ------------------------------
LICDIR="$HOME/.matlab/R2026a_licenses"; mkdir -p "$LICDIR"
if ! ls "$LICDIR"/*.lic >/dev/null 2>&1 && [ -z "${MLM_LICENSE_FILE:-}" ]; then
  log "NOTE: no MATLAB license detected. Provision it out of band before running instances:"
  log "        file license:    copy the .lic into $LICDIR/"
  log "        network license: export MLM_LICENSE_FILE=port@host"
  log "        online license:  the MathWorks AMI signs in on first GUI launch (flaky across"
  log "                         stop/start -- a clean restart re-checks-out the token)"
fi

# ---- 4. NNV deps + ONNX/PyTorch converters + python + matlab.engine (maintained installer) -
INSTALL_TOOL="$HOME/nnv/code/nnv/examples/Submission/VNN_COMP2026/install_tool.sh"
if [ -f "$INSTALL_TOOL" ]; then
  log "running install_tool.sh v1 (mpm converters, NNV install, onnx/onnxruntime, matlab.engine)"
  sudo bash "$INSTALL_TOOL" v1 || log "WARN: install_tool.sh returned non-zero -- inspect output above"
else
  log "ERROR: install_tool.sh not found at $INSTALL_TOOL (wrong NNV_BRANCH?)"
fi

# ---- 5. smoke check -----------------------------------------------------------------------
log "smoke: NNV version + ONNX importer availability"
matlab -batch "cd('$HOME/nnv/code/nnv'); startup_nnv; fprintf('NNV %s\n', NNVVERSION()); fprintf('importNetworkFromONNX present: %d\n', exist('importNetworkFromONNX'));" 2>&1 | tail -6 || log "WARN: smoke failed"
log "DONE. nnv=~/nnv ($NNV_BRANCH); benchmarks 2025+2026; results 2025; status repo cloned."
log "Next: bash ~/nnv/code/nnv/examples/Submission/VNN_COMP2026/aws/overnight_sweep.sh"
