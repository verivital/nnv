#!/bin/bash

# Installation script used for VNN-COMP. Ubuntu, MATLAB R2022b (https://github.com/verivital/nnv)

TOOL_NAME="nnv"
VERSION_STRING="v1"

# check arguments
if [ "$1" != ${VERSION_STRING} ]; then
	echo "Expected first argument (version string) '$VERSION_STRING', got '$1'"
	exit 1
fi

echo "Installing $TOOL_NAME dependencies"

sudo apt install -y python3-pip
pip install numpy onnxruntime onnx scipy matlabengine

# remove paths from any prior installation
matlab -nodisplay -r "rmpath(genpath('/home/ubuntu/work/nnv/')); savepath; quit"

matlab -nodisplay -r "cd '/home/ubuntu/work/nnv/code/nnv/'; install; addpath(genpath('/home/ubuntu/work/nnv/')); savepath; quit"

# get paths to ensure support packages are set properly
matlab -nodisplay -r "matlabshared.supportpkg.getInstalled; matlabshared.supportpkg.getSupportPackageRoot; matlabroot; ver; quit"
