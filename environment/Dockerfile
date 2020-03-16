# hash:sha256:f622d7353e92f49ac2e11234b75ba3e52d96df0222aa7e2bf36ad6cd44b3be51
FROM registry.codeocean.com/codeocean/matlab:2019a-ubuntu18.04

# codeocean does not seem to support multi-stage builds, so python set up is below versus using this image
# hash:sha256:1f62953520a0e0fb91a3b47349bfbecf2b530ada4a23f795534f50d6d9db554b
#FROM registry.codeocean.com/codeocean/miniconda3:4.7.10-cuda10.1-cudnn7-ubuntu18.04

# older matlab version (2017b)
# hash:sha256:3e64c2baee87bdfdf6e7bac9b76b188098fce8e831ecac85fd33a589aa3c1419
#FROM registry.codeocean.com/codeocean/matlab:2017b-ubuntu16.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

#RUN apt-add-repository universe

RUN apt-get update \
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

RUN apt-get install -y texlive-full
RUN apt-get install -y nano

RUN pip install numpy
RUN pip3 install numpy
RUN pip install pybind11
RUN pip3 install pybind11

RUN pkg-config --cflags python

# set up comparison scripts (for ReluVal, Reluplex, Marabou)

# this directory is used as a flag to detect if on codeocean or not
RUN mkdir -p /codeocean-true

RUN mkdir -p /deps

RUN ls /codeocean-true

WORKDIR /deps

RUN git clone https://github.com/verivital/nnvmt.git

# these ls commands are just to force a new build so that the cloned repo is latest and not the cached one in the prior docker container image
RUN ls
RUN ls

RUN git clone https://github.com/verivital/run_nnv_comparison.git
WORKDIR /deps/run_nnv_comparison 
RUN ls
RUN chmod +x setup_tools.sh

# setting up all the tools will take several minutes at least, so comment out if comparison is not needed
RUN ./setup_tools.sh

WORKDIR /

COPY postInstall /
RUN /postInstall

# set up matlab support packages (for vgg)
WORKDIR /MATLAB/bin/glnxa64
RUN ls
#RUN ./SupportSoftwareInstaller -downloadfolder /data/matlab_sp -inputfile /data/matlab_sp/ssi_input.txt
# previous gives error:
# 
#root@d5a9c3095b6b:/MATLAB/bin/glnxa64# ./SupportSoftwareInstaller
#terminate called after throwing an instance of '(anonymous namespace)::DisplayError'
#  what():  No display available.
#Aborted (core dumped)
