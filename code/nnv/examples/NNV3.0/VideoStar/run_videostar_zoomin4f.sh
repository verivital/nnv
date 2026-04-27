#!/bin/bash
# VideoStar ZoomIn-4f Verification Script
#
# Thin wrapper around run_zoomin_4f.m. The .m file is the canonical
# reference implementation; this wrapper is convenient for shell-only
# automation and CI.
#
# Usage:
#   ./run_videostar_zoomin4f.sh            # from this directory
#
# Equivalent direct invocation:
#   matlab -batch "run_zoomin_4f"
#
# A Python wrapper (run_zoomin_4f.py) is also provided.

set -euo pipefail
cd "$(dirname "$0")"

if [ ! -f "run_zoomin_4f.m" ]; then
    echo "Error: run_zoomin_4f.m not found in $(pwd)" >&2
    exit 1
fi

echo "=============================================="
echo "VideoStar ZoomIn-4f Verification"
echo "=============================================="

matlab -nodisplay -batch "run_zoomin_4f"

echo ""
echo "=============================================="
echo "Done. Results under VideoStar/results/<timestamp>/"
echo "=============================================="
