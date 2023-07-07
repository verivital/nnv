#!/bin/bash

# example run_instance.sh script for VNNCOMP 2023 for nnv (https://github.com/verivital/nnv.git)
# four arguments, first is "v1", second is a benchmark category identifier string such as "acasxu", third is a path to the .onnx file and fourth is a path to .vnnlib file
# modified: Samuel Sasaki, June 28th 2023

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
#DIR=$(dirname $(dirname $(realpath $0)))
#export PYTHONPATH="$PYTHONPATH:$DIR/src"

echo $PYTHONPATH
# TIMEOUT=$5

# echo "Running $TOOL_NAME on benchmark instance in category '$CATEGORY' with onnx file '$ONNX_FILE', vnnlib file '$VNNLIB_FILE', and timeout $TIMEOUT"

python3 /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2023/execute.py 'run_instance' $CATEGORY $ONNX_FILE $VNNLIB_FILE $TIMEOUT $RESULTS_FILE
# sudo python3 execute.py 'run_instance' "$CATEGORY" "$ONNX_FILE" "$VNNLIB_FILE" "$TIMEOUT" "$RESULTS_FILE"

echo ""
