#!/bin/bash

# example prepare_instance.sh script for VNNCOMP 2023 for nnv (https://github.com/verivital/nnv.git)
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

echo "Preparing $TOOL_NAME for benchmark instance in category '$CATEGORY' with onnx file '$ONNX_FILE' and vnnlib file '$VNNLIB_FILE'"

# kill any zombie processes
killall -q python3
# killall -q python
killall -q matlab

#WAIT_FOR_CONNECTION_TO_CLOSE='import matlab.engine\nimport time\nwhile matlab.engine.find_matlab(): time.sleep(1)'
# python3 -c "exec('$WAIT_FOR_CONNECTION_TO_CLOSE')"
# python -c "exec('$WAIT_FOR_CONNECTION_TO_CLOSE')"

# start the matlab engine in background and keep the connection open
# python3 /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2023/execute.py 'prepare_instance' $CATEGORY $ONNX_FILE $VNNLIB_FILE &
# python execute.py 'prepare_instance' $CATEGORY $ONNX_FILE $VNNLIB_FILE &

WAIT_FOR_CONNECTION_TO_OPEN='import matlab.engine\nimport time\nwhile not matlab.engine.find_matlab(): time.sleep(1) \nprint(eng)'
# python3 -c "exec('$WAIT_FOR_CONNECTION_TO_OPEN')"
# python -c "exec('$WAIT_FOR_CONNECTION_TO_OPEN')"

exit 0
