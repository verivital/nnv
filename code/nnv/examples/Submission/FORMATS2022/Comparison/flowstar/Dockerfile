FROM ubuntu:16.04

ENV FLOWSTAR_VERSION 2.1.0
ENV FLOWSTAR_FILE_SHA512SUM 'd5243f3bbcdd6bffcaf2f1ae8559278f62567877021981e4443cd90fbf2918e0acb317a2d27724bc81d3a0e38ad7f7d48c59d680be1dd5345e80d2234dd3fe3b'

RUN mkdir -p /flowstar
WORKDIR /flowstar
RUN apt-get update && apt-get install make
RUN apt-get install -qy curl flex bison libgmp-dev libmpfr-dev libgsl-dev gnuplot
RUN apt-get install -y make g++ wget openjdk-8-jdk unzip python3 \
    && apt-get install -y libgmp3-dev libmpfr-dev libmpfr-doc libmpfr4 libmpfr4-dbg gsl-bin libgsl0-dev bison flex gnuplot-x11 libglpk-dev \
    && rm -rf /var/lib/apt/lists/*
RUN curl -fL https://github.com/stanleybak/hybrid_tools/raw/master/flowstar-${FLOWSTAR_VERSION}.tar.gz > flowstar.tar.gz
# print and check hash
RUN sha512sum flowstar.tar.gz | tee flowstar.sha512sum && grep -q "${FLOWSTAR_FILE_SHA512SUM}" flowstar.sha512sum
RUN tar xzf flowstar.tar.gz
WORKDIR /flowstar/flowstar-${FLOWSTAR_VERSION}/
RUN make
# ENV PATH=$PATH:/flowstar-${FLOWSTAR_VERSION}/
# ADD Benchmarks /flowstar/flowstar-${FLOWSTAR_VERSION}//Benchmarks

# copy current directory into container
COPY . /flowstar/flowstar-2.1.0/

RUN chmod +x run_all.sh
