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

# ONNX importer is useless with command line...
# curl --retry 100 --retry-connrefused -L -O https://www.dropbox.com/s/4p17xm4tlm8r9gs/sppFile.zip
# sleep 60
#  sudo mkdir -p /usr/local/MATLAB/R2022b/SupportPackages
# sudo unzip sppFile.zip -d /usr/local/MATLAB/R2022b/SupportPackages

# LOOKS LIKE WE STILL NEED ONNX BECAUSE SOME LAYERS ARE IN THERE LIKE nnet.onnx.layer.ElementwiseAffineLayer
# Does adding those functions to the path without actually installing them work??
curl --retry 100 --retry-connrefused -L -O https://www.dropbox.com/s/4p17xm4tlm8r9gs/sppFile.zip
sleep 60
sudo unzip sppFile.zip -d cd /home/ubuntu/toolkit/code/nnv
sudo rm sppFile.zip

cd /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2023/
curl --retry 100 --retry-connrefused -L -O https://www.dropbox.com/scl/fi/7wr6dy1qe7k8c670ieunb/networks2023.zip?rlkey=eo40pf96qtvdhz11civ1u6jt1&dl=0

sleep 60

unzip  *.zip*

ip link show # get mac address (for licensing)

echo $USER # get usernme (for licensing)

mkdir ~/.matlab/R2022b_licenses 