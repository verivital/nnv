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
DIR=$(dirname $(dirname $(dirname $(dirname $(dirname $(realpath $0))))))

mkdir -p $DIR/deps
cd deps 

sudo apt install python3-pip
pip3 install numpy
pip3 install numpy
pip3 install pybind11

apt-get install -y libprotobuf-dev protobuf-compiler

apt-get install -y psmisc && # for killall, used in prepare_instance.sh script
pip3 install -r "$DIR/dep_requirements.txt"

mkdir -p $DIR/tbxmanager && cd tbxmanager
> install.m

matlab -nodisplay -r "install; addpath(genpath('../../deps')); savepath"

matlab -nodisplay -r "addpath(genpath('$DIR/deps')); addpath(genpath('$DIR./code')); savepath;"


