#!/bin/bash

# Installation script used for VNN-COMP. Ubuntu, MATLAB R2024b (https://github.com/verivital/nnv)

TOOL_NAME="nnv"
VERSION_STRING="v1"

# check arguments
if [ "$1" != ${VERSION_STRING} ]; then
	echo "Expected first argument (version string) '$VERSION_STRING', got '$1'"
	exit 1
fi

echo "Installing $TOOL_NAME dependencies"

#ip link show # get mac address (for licensing)

#echo $USER # get username (for licensing)

#mkdir ~/.matlab/R2024b_licenses 

# INSTALL MPM TO INSTALL ADDITIONAL MATLAB PACKAGES

apt install wget

wget https://www.mathworks.com/mpm/glnxa64/mpm

chmod +x mpm

#./mpm install --release=R2024b --products Deep_Learning_Toolbox_Converter_for_ONNX_Model_Format --destination /home/ubuntu/Documents/MATLAB/SupportPackages/R2024b

./mpm install --release=R2024b --products Deep_Learning_Toolbox_Converter_for_ONNX_Model_Format --destination /usr/local/matlab

# ADD STEPS TO INSTALL GUROBI (optional)

#cd ~/

#wget https://packages.gurobi.com/11.0/gurobi11.0.2_linux64.tar.gz

#tar xvfz gurobi11.0.2_linux64.tar.gz

