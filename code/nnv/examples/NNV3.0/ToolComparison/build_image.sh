#!/bin/bash
# build_image.sh — build nnv3.0:r2025b with the Vanderbilt MATLAB license.
#
# Idempotent thanks to Docker's layer cache. First build is 30-60 min;
# subsequent rebuilds (when host code changes) are minutes.
#
# After build, run_all.sh expects the image tag nnv3.0:r2025b.

set -euo pipefail

IMAGE_TAG="nnv3.0:r2025b"
LICENSE_SERVER="${LICENSE_SERVER:-27009@licenseserver.it.vanderbilt.edu}"
HOST_REPO="$(cd "$(dirname "$0")/../../../../.." && pwd)"
DOCKERFILE="code/nnv/examples/NNV3.0/Dockerfile"

echo "=== Docker build ==="
echo "  tag:        $IMAGE_TAG"
echo "  Dockerfile: $DOCKERFILE"
echo "  license:    $LICENSE_SERVER"
echo "  context:    $HOST_REPO"
echo

cd "$HOST_REPO"
docker build \
    -t "$IMAGE_TAG" \
    -f "$DOCKERFILE" \
    --build-arg MATLAB_RELEASE=R2025b \
    --build-arg LICENSE_SERVER="$LICENSE_SERVER" \
    .

echo
echo "=== Image built. Inspect: ==="
docker image inspect "$IMAGE_TAG" --format '{{.RepoTags}} {{.Size}} bytes'
