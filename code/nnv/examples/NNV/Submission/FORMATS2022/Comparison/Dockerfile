FROM ubuntu:20.04

RUN apt-get update && \
    apt-get install -y python3.8 python3.8 python3-pip && \
    apt-get clean

# Install everything with no interactions   
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago

RUN apt-get install -qy curl unzip

RUN curl -fL https://github.com/mldiego/GoTube/archive/refs/heads/main.zip > gotube.zip

RUN unzip gotube.zip

WORKDIR GoTube-main

RUN chmod +x paper_submission_compare.sh

RUN pip3 install -r requirements.txt

RUN chmod +x paper_submission_short.sh

# Clear previous results from GitHub
RUN mv /GoTube-main/logged /GoTube-main/orig_logged
RUN mv /GoTube-main/saved_outputs /GoTube-main/orig_saved_outputs

# Installing Flow* now
WORKDIR ..

# FROM ubuntu:16.04

ENV FLOWSTAR_VERSION 2.1.0
ENV FLOWSTAR_FILE_SHA512SUM 'd5243f3bbcdd6bffcaf2f1ae8559278f62567877021981e4443cd90fbf2918e0acb317a2d27724bc81d3a0e38ad7f7d48c59d680be1dd5345e80d2234dd3fe3b'

RUN mkdir -p /flowstar
WORKDIR /flowstar
RUN apt-get update && apt-get install make
#RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt-get install -qy curl flex bison libgmp-dev libmpfr-dev libgsl-dev gnuplot
RUN apt-get install -y make g++ wget openjdk-8-jdk unzip python3 \
    && apt-get install -y libgmp3-dev libmpfr-dev libmpfr-doc gsl-bin libgsl0-dev bison flex gnuplot-x11 libglpk-dev \
    && rm -rf /var/lib/apt/lists/*
RUN curl -fL https://github.com/stanleybak/hybrid_tools/raw/master/flowstar-${FLOWSTAR_VERSION}.tar.gz > flowstar.tar.gz
# print and check hash
RUN sha512sum flowstar.tar.gz | tee flowstar.sha512sum && grep -q "${FLOWSTAR_FILE_SHA512SUM}" flowstar.sha512sum
RUN tar xzf flowstar.tar.gz
WORKDIR /flowstar/flowstar-${FLOWSTAR_VERSION}/
RUN make

# copy current directory into container
COPY flowstar/ /flowstar/flowstar-2.1.0/

RUN chmod +x run_all.sh
RUN chmod +x run_short.sh

RUN mv /flowstar/flowstar-2.1.0/results /flowstar/flowstar-2.1.0/orig_results

WORKDIR ../..

# Installing Julia and Juliareach
RUN apt-get update && apt-get -qy install make gcc hdf5-tools libsz2 libz-dev libqt5widgets5 texlive-latex-base dvipng

RUN apt -y install build-essential

RUN apt-get install wget

RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.6/julia-1.6.0-linux-x86_64.tar.gz
RUN tar xvfz julia-1.6.0-linux-x86_64.tar.gz
RUN ln -s /home/ubuntu/julia-1.6.0/bin/julia /usr/local/bin/julia

# set working directory
WORKDIR /juliareach

# copy juliareach directory into container
COPY juliareach/ /juliareach

RUN mv /juliareach/results /juliareach/orig_results 
# add julia to the path
ENV PATH="$PATH:/julia-1.6.0/bin"
WORKDIR ..
# Copy run_tools files
COPY run_tools /run_tools
COPY run_tools_short /run_tools_short


