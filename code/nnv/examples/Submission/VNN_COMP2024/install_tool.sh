#!/bin/bash

# Installation script used for VNN-COMP. Ubuntu, MATLAB R2024a (https://github.com/verivital/nnv)

TOOL_NAME="nnv"
VERSION_STRING="v1"

# check arguments
if [ "$1" != ${VERSION_STRING} ]; then
	echo "Expected first argument (version string) '$VERSION_STRING', got '$1'"
	exit 1
fi

echo "Installing $TOOL_NAME dependencies"

# Install support packages (ONNX importer)
curl --retry 100 --retry-connrefused -L -O https://www.dropbox.com/s/4p17xm4tlm8r9gs/sppFile.zip # need to do this step for R2024a
sleep 60
sudo unzip sppFile.zip -d /home/ubuntu/toolkit/code/nnv
sudo rm sppFile.zip

# This is called networks2023, but it contains the MATLAB networks for most of 2022 and 2023 (Do we need this step for 2024?)
cd /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2023/
curl --retry 100 --retry-connrefused -L -O https://www.dropbox.com/scl/fi/chb0fotern4r5sycz57r9/networks2023.zip?rlkey=1cwvh73b5zgtnb67xdrvlg8am&dl=0

sleep 60

unzip  *.zip*

ip link show # get mac address (for licensing)

echo $USER # get usernme (for licensing)

mkdir ~/.matlab/R2024a_licenses 

# ADD STEPS TO INSTALL GUROBI