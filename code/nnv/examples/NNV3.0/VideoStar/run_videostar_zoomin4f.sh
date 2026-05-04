#!/bin/bash
# VideoStar ZoomIn-4f Verification Script
#
# This script runs a subset of the ZoomIn-4f verification experiments
# for VideoStar. It is designed to be run inside the NNV Docker container.
#
# Usage (inside Docker container):
#   cd /home/matlab/nnv/code/nnv/examples/NNV3.0/VideoStar
#   ./run_videostar_zoomin4f.sh
#
# Or with MATLAB headless:
#   cd /home/matlab/nnv/code/nnv/examples/NNV3.0/VideoStar
#   matlab -nodisplay -r "run('run_zoomin_4f.m'); exit()"
#
# Or with Python:
#   cd /home/matlab/nnv/code/nnv/examples/NNV3.0/VideoStar
#   python run_zoomin_4f.py --algorithm relax --num-samples 10

echo "=============================================="
echo "VideoStar ZoomIn-4f Verification"
echo "=============================================="
echo ""

# Check if we're in the right directory
if [ ! -f "run_zoomin_4f.m" ]; then
    echo "Error: run_zoomin_4f.m not found."
    echo "Please run this script from the VideoStar directory:"
    echo "  cd /home/matlab/nnv/code/nnv/examples/NNV3.0/VideoStar"
    exit 1
fi

echo "Running ZoomIn-4f verification with MATLAB..."
echo ""

# Run the MATLAB script
matlab -nodisplay -r "run('run_zoomin_4f.m'); exit()"

echo ""
echo "=============================================="
echo "VideoStar ZoomIn-4f Verification Complete"
echo "=============================================="
echo ""
echo "Results saved to: /tmp/results/VideoStar/ZoomIn/4/"
echo ""
echo "To copy results from the container:"
echo "  docker cp <container_id>:/tmp/results/VideoStar ./videostar_results"
