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

curl --retry 100 --retry-connrefused -L -O https://www.dropbox.com/s/4p17xm4tlm8r9gs/sppFile.zip

sleep 60

sudo mkdir -p /usr/local/MATLAB/R2022b/SupportPackages
sudo unzip sppFile.zip -d /usr/local/MATLAB/R2022b/SupportPackages

ip link show # get mac address (for licensing)

echo $USER # get usernme (for licensing)

mkdir ~/.matlab/R2022b_licenses 