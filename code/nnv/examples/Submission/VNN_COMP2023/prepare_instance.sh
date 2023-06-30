#!/bin/bash
# example prepare_instance.sh script for VNNCOMP2023 for nnv (https://github.com/verivital/nnv.git) 

# four arguments, first is "v1", second is a benchmark category identifier string such as "acasxu", third is path to the .onnx file and fourth is path to .vnnlib file

TOOL_NAME="nnv"
VERSION_STRING="v1"

# # check arguments
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
killall -q matlab


exit 0
