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

apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        gcc \
        libboost-all-dev \
        libopenblas-base \
        make \
        pkg-config \
        python3 \
        python3-dev \
        wget \
        python-pip \
        python3-pip \
        software-properties-common

apt-get install -y texlive-full
apt-get install -y nano

pip install numpy
pip3 install numpy
pip install pybind11
pip3 install pybind11

pkg-config --cflags python

apt-get install -y libprotobuf-dev protobuf-compiler

apt-get install -y psmisc && # for killall, used in prepare_instance.sh script
pip3 install -r "$DIR/dep_requirements.txt"

curl -LO https://tumcps.github.io/CORA/data/CORA_2018.zip
unzip CORA_2018.zip && rm CORA_2018.zip

mkdir -p $DIR/tbxmanager && cd tbxmanager
> install.m

matlab -nodisplay -r "install; addpath(genpath('../../deps')); savepath"

matlab -nodisplay -r "addpath(genpath('$DIR/deps')); addpath(genpath('$DIR./code')); savepath;"

cd /MATLAB/extern/engines/python
python3 setup.py install

