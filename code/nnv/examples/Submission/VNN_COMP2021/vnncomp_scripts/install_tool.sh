#!/bin/bash
# install_tool.sh script for VNNCOMP2021 for nnv (https://github.com/verivital/nnv.git) 

# modified: Neelanjana , June 25, 2021

TOOL_NAME="nnv"
VERSION_STRING="v1"

# check arguments
if [ "$1" != ${VERSION_STRING} ]; then
	echo "Expected first argument (version string) '$VERSION_STRING', got '$1'"
	exit 1
fi

echo "Installing $TOOL_NAME dependencies"

sudo apt install -y python3-pip
pip install numpy onnxruntime onnx scipy

# remove paths from any prior installation
matlab -nodisplay -r "rmpath(genpath('/home/ubuntu/work/nnv/')); savepath; quit"

matlab -nodisplay -r "cd '/home/ubuntu/work/nnv/code/nnv/'; install; addpath(genpath('/home/ubuntu/work/nnv/')); savepath; quit"



