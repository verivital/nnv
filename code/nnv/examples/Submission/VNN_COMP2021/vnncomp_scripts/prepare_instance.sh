#!/bin/bash
# example prepare_instance.sh script for VNNCOMP2021 for nnv (https://github.com/verivital/nnv.git) 

# four arguments, first is "v1", second is a benchmark category identifier string such as "acasxu", third is path to the .onnx file and fourth is path to .vnnlib file

# modified: Neelanjana , June 25, 2021

TOOL_NAME= nnv
VERSION_STRING=v1

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
pgrep -f matlab | xargs kill -9

# script returns a 0 exit code if successful. If you want to skip a benchmark category you can return non-zero.
if [ "$2" == 'cifar10_resnet' ] || [ "$2" == 'nn4sys' ] || [ "$2" == 'marabou-cifar10' ] || [ "$2" == 'test' ] || || [ "$2" == 'mnistfc' ]; then
	echo "NNV is not participating in category '$CATEGORY' " 
	exit 2
fi

python3 ../get_specs.py "$ONNX_FILE" "$VNNLIB_FILE"

matlab -nodisplay -r "addpath(genpath('../../../../../../code')); preProcessing('$ONNX_FILE','$VNNLIB_FILE');exit;" 

# kill any zombie processes
killall -q python3
pgrep -f matlab | xargs kill -9

exit 0
