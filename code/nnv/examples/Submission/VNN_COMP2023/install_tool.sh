#!/bin/bash

# Installation script used for VNN-COMP. Ubuntu, MATLAB R2022b (https://github.com/verivital/nnv)

TOOL_NAME="nnv"
VERSION_STRING="v1"

# check arguments
if [ "$1" != ${VERSION_STRING} ]; then
	echo "Expected first argument (version string) '$VERSION_STRING', got '$1'"
	exit 1
fi

# sudo killall apt apt-get

echo "Installing $TOOL_NAME dependencies"

sudo apt install -y python3-pip 
# pip install numpy matlabengine
#echo "Finding MATLAB (directory)"
find / -type d -iname "MATLAB" # is matlab installed? where?
#echo "Finding MATLAB (anything)"
find / -type f -iname "MATLAB" # is matlab installed? where?

find / -type f -iname "license"
find / -type d -iname "license"
find / -type f -name "*.lic"

# ifconfig

# remove paths from any prior installation
matlab -nodisplay -r "rmpath(genpath('/home/ubuntu/work/nnv/')); savepath; quit"

matlab -r "matlabshared.supportpkg.getInstalled; matlabshared.supportpkg.getSupportPackageRoot; matlabroot; ver; quit"

#matlab -nodisplay -r "cd '/home/ubuntu/work/nnv/code/nnv/'; install; addpath(genpath('/home/ubuntu/work/nnv/')); savepath; quit"

# get paths to ensure support packages are set properly
#matlab -nodisplay -r "matlabshared.supportpkg.getInstalled; matlabshared.supportpkg.getSupportPackageRoot; matlabroot; ver; quit"
