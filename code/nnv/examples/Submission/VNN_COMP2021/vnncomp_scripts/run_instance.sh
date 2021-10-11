#!/bin/bash
# example run_benchmark.sh script for VNNCOMP2021 for nnv (https://github.com/verivital/nnv.git) 

# six arguments, first is "v1", second is a benchmark category itentifier string such as "acasxu", third is path to the .onnx file, fourth is path to .vnnlib file, fifth is a path to the results file, and sixth is a timeout in seconds.

# modified: Neelanjana , June 25, 2021 

TOOL_NAME="nnv"
VERSION_STRING="v1"

# check arguments
if [ "$1" != ${VERSION_STRING} ]; then
	echo "Expected first argument (version string) '$VERSION_STRING', got '$1'"
	exit 1
fi

CATEGORY=$2
ONNX_FILE=$3
VNNLIB_FILE=$4
RESULTS_FILE=$5
TIMEOUT=$6

echo "Running $TOOL_NAME on benchmark instance in category '$CATEGORY' with onnx file '$ONNX_FILE', vnnlib file '$VNNLIB_FILE', results file $RESULTS_FILE, and timeout $TIMEOUT"

# setup environment variable for tool (doing it earlier won't be persistent with docker)"
DIR=$(dirname $(dirname $(realpath $0)))
export PYTHONPATH="$PYTHONPATH:$DIR/src"

#export OPENBLAS_NUM_THREADS=1
#export OMP_NUM_THREADS=1

# run the tool to produce the results file
python3 /home/ubuntu/work/nnv/code/nnv/examples/Submission/VNN_COMP2021/execute_runs.py "$ONNX_FILE" "$VNNLIB_FILE" "$TIMEOUT" "$RESULTS_FILE" "CATEGORY"
