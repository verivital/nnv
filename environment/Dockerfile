# hash:sha256:f0637c44f33f41531b237d9548759f2f321a1c02db722190f4f0d27932dab985
FROM registry.codeocean.com/codeocean/matlab:2022b-ubuntu20.04

# old 2021b version
# hash:sha256:5f616c6f66b2c7c0f5e3b8d5b21214210b9d8d1dbe6b703011b15b8d65208c91
#FROM registry.codeocean.com/codeocean/matlab:2021b-ubuntu20.04

# old 2020b version
# hash:sha256:eb133fa6ef390c8a024c045858d448a3c7a34f140ef89a10c693182c8e4e31c7
#FROM registry.codeocean.com/codeocean/matlab:2020b-ubuntu20.04

# old 2019a version
## hash:sha256:f622d7353e92f49ac2e11234b75ba3e52d96df0222aa7e2bf36ad6cd44b3be51
#FROM registry.codeocean.com/codeocean/matlab:2019a-ubuntu18.04

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
        python3-pip \
        software-properties-common

#RUN apt-get install -y texlive-full # used for regen some plots, but very long install, avoid if possible
RUN apt-get install -y nano

RUN pip install numpy
RUN pip3 install numpy
RUN pip install pybind11
RUN pip3 install pybind11

#RUN pkg-config --cflags python

# set up comparison scripts (for ReluVal, Reluplex, Marabou)

# this directory is used as a flag to detect if on codeocean or not
RUN mkdir -p /codeocean-true

RUN mkdir -p /deps

RUN ls /codeocean-true


# start of formats setup for comparisons
# notes: may need to be careful with python versions, is the 3.8 needed?
# Comment out unless reproducing formats2022 results
# RUN apt-get update && \
#    apt-get install -y python3.8 python3.8 python3-pip && \
#    apt-get clean

# Install everything with no interactions   
#ARG DEBIAN_FRONTEND=noninteractive
#ENV TZ=America/Chicago

#RUN apt-get install -qy curl unzip

#RUN curl -fL https://github.com/mldiego/GoTube/archive/refs/heads/main.zip > gotube.zip

#RUN unzip gotube.zip

#WORKDIR GoTube-main

#RUN chmod +x paper_submission_compare.sh

#RUN pip install -r requirements.txt
#RUN pip3 install -r requirements.txt

#RUN chmod +x paper_submission_short.sh

# Clear previous results from GitHub
#RUN mv /GoTube-main/logged /GoTube-main/orig_logged
#RUN mv /GoTube-main/saved_outputs /GoTube-main/orig_saved_outputs

# Installing Flow* now
#WORKDIR ..

# FROM ubuntu:16.04

#ENV FLOWSTAR_VERSION 2.1.0
#ENV FLOWSTAR_FILE_SHA512SUM 'd5243f3bbcdd6bffcaf2f1ae8559278f62567877021981e4443cd90fbf2918e0acb317a2d27724bc81d3a0e38ad7f7d48c59d680be1dd5345e80d2234dd3fe3b'

#RUN mkdir -p /flowstar
#WORKDIR /flowstar
#RUN apt-get update && apt-get install make
#RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
#RUN apt-get install -qy curl flex bison libgmp-dev libmpfr-dev libgsl-dev gnuplot
#RUN apt-get install -y make g++ wget openjdk-8-jdk unzip python3 \
#    && apt-get install -y libgmp3-dev libmpfr-dev libmpfr-doc gsl-bin libgsl0-dev bison flex gnuplot-x11 libglpk-dev \
#    && rm -rf /var/lib/apt/lists/*
#RUN curl -fL https://github.com/stanleybak/hybrid_tools/raw/master/flowstar-${FLOWSTAR_VERSION}.tar.gz > flowstar.tar.gz
# print and check hash
#RUN sha512sum flowstar.tar.gz | tee flowstar.sha512sum && grep -q "${FLOWSTAR_FILE_SHA512SUM}" flowstar.sha512sum
#RUN tar xzf flowstar.tar.gz
#WORKDIR /flowstar/flowstar-${FLOWSTAR_VERSION}/
#RUN make

# copy current directory into container
#RUN ls
#RUN ls /flowstar
# this copy is from the git repository folder (on host) to the container; shouldn't be necessary as we're already in the host with this mounted elsewhere already, and it is causing an error; so drop it
#COPY flowstar/ /flowstar/flowstar-2.1.0/

#RUN chmod +x run_all.sh
#RUN chmod +x run_short.sh

#RUN mv /flowstar/flowstar-2.1.0/results /flowstar/flowstar-2.1.0/orig_results

#WORKDIR ../..

# Installing Julia and Juliareach
#RUN apt-get update && apt-get -qy install make gcc hdf5-tools libsz2 libz-dev libqt5widgets5 texlive-latex-base dvipng

#RUN apt -y install build-essential

#RUN apt-get install wget

#RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.0-linux-x86_64.tar.gz
#RUN tar xvfz julia-1.6.0-linux-x86_64.tar.gz
#RUN ln -s /home/ubuntu/julia-1.6.0/bin/julia /usr/local/bin/julia

# set working directory
#WORKDIR /juliareach

# copy juliareach directory into container
# same as with flowstar: already available in container, so don't do this
#COPY juliareach/ /juliareach

#RUN ls
#RUN ls /juliareach

#RUN mv /juliareach/results /juliareach/orig_results 
# add julia to the path
#ENV PATH="$PATH:/julia-1.6.0/bin"
#WORKDIR ..
# Copy run_tools files
#COPY run_tools /run_tools
#COPY run_tools_short /run_tools_short
# end of formats setup for comparison



WORKDIR /deps

RUN git clone https://github.com/verivital/nnvmt.git

# handled via postInstall zip release sometimes, but we need some bug fixes from latest release that we made for FORMATS'22 package
RUN git clone https://github.com/verivital/CORA.git

# these ls commands are just to force a new build so that the cloned repo is latest and not the cached one in the prior docker container image

RUN git clone https://github.com/verivital/run_nnv_comparison.git
WORKDIR /deps/run_nnv_comparison 
RUN ls
RUN chmod +x setup_tools.sh

# setting up all the tools will take several minutes at least, so comment out if comparison is not needed
#RUN ./setup_tools.sh

WORKDIR /deps/nnvmt
RUN ln -s /code/nnv/nnmt /deps/nnvmt 

# nnmt needs anaconda
# #https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
#RUN wget https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh
#RUN bash Anaconda3-2020.02-Linux-x86_64.sh -b
#RUN rm Anaconda3-2020.02-Linux-x86_64.sh

# Set path to conda
#ENV PATH /root/anaconda3/bin:$PATH
#ENV PATH /home/ubuntu/anaconda3/bin:$PATH

RUN ls /
RUN ls /home/
RUN ls /usr
RUN ls /root

# Updating Anaconda packages
#RUN conda update conda
#RUN conda update anaconda
#RUN conda update --all


#WORKDIR /deps/nnvmt
#RUN conda create -n nnmt python=3.6
#RUN conda init bash
#RUN source ~/anaconda3/etc/profile.d/conda.sh
#RUN conda activate nnmt

RUN apt-get install -y libprotobuf-dev protobuf-compiler

#RUN pip install --user -r requirements.txt

# there was an error related to version differences: 
#RUN pip install --upgrade onnx

WORKDIR /

COPY postInstall /
RUN /postInstall

# set up matlab support packages (for vgg)
#WORKDIR /MATLAB/bin/glnxa64
#RUN ls
#RUN ./SupportSoftwareInstaller -downloadfolder /data/matlab_sp -inputfile /data/matlab_sp/ssi_input.txt
# previous gives error:
# 
#root@d5a9c3095b6b:/MATLAB/bin/glnxa64# ./SupportSoftwareInstaller
#terminate called after throwing an instance of '(anonymous namespace)::DisplayError'
#  what():  No display available.
#Aborted (core dumped)
