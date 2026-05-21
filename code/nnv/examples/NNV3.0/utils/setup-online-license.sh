#!/usr/bin/env bash
# setup-online-license.sh -- one-time NNV + AIVL installation for the
# nnv3.0-online image.
#
# Run this INSIDE the running nnv3.0-online container, ONCE, after you
# have signed in to MATLAB via the browser at http://localhost:8888.
# It runs the two MATLAB invocations that the standard (network-licence)
# Dockerfile does at build time -- here they have to wait until a
# licence is available, which is after sign-in.
#
# Subsequent `docker run` invocations of nnv3.0-online with the same
# named volumes mounted (`-v nnv3-matlab-prefs:/home/matlab/.matlab
# -v nnv3-matlab-mw:/home/matlab/.MathWorks`) reuse this setup; you
# do not need to re-run unless the volumes are removed or the marker
# file below is deleted.
#
# Usage from the host:
#   docker exec -it <container-name> \
#     bash /home/matlab/nnv/code/nnv/examples/NNV3.0/utils/setup-online-license.sh

set -euo pipefail

MARKER=/home/matlab/.matlab/.nnv-online-setup-done
NNV_ROOT=/home/matlab/nnv
VENV_PY=${NNV_ROOT}/.venv/bin/python

if [ -f "$MARKER" ]; then
    echo "[setup] NNV + AIVL already installed in this volume."
    echo "[setup] Marker: $MARKER"
    echo "[setup] Delete the marker to force re-run."
    exit 0
fi

echo "[setup] (1/2) Installing NNV paths and verifying setup..."
matlab -nodisplay -batch "\
    pyenv('Version', '${VENV_PY}'); \
    cd('${NNV_ROOT}/code/nnv'); \
    try, install; catch ME, fprintf(2, '[install warning] %s\\n', ME.message); end; \
    savepath; \
    cd('${NNV_ROOT}/code/nnv'); check_nnv_setup(); \
    "

echo "[setup] (2/2) Extracting AIVL Support Package (skipped if tarball missing)..."
matlab -nodisplay -batch "\
    cd('${NNV_ROOT}/code/nnv/examples/NNV3.0/ToolComparison/utils'); \
    try, run('toolbox_install.m'); catch ME, fprintf(2, '[toolbox_install warning] %s\\n', ME.message); end; \
    "

mkdir -p "$(dirname "$MARKER")"
touch "$MARKER"

echo "[setup] Done."
echo "[setup] Future runs of nnv3.0-online with the same named volumes"
echo "[setup] mounted will reuse this NNV install and the cached MATLAB"
echo "[setup] activation. Exit the interactive container (Ctrl-C) and run"
echo "[setup] experiments headlessly per the Dockerfile header comment."
