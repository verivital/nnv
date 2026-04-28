# NNV3.0 — repeatability-package Dockerfile.
#
# Builds a MATLAB R2024b image with the NNV repo and the NNV3.0
# example assets baked in. The build is **license-agnostic** — MATLAB
# is installed via mpm but never launched, so no license check happens
# at build time. License is supplied at `docker run` time. See the
# top-level README's "License setup" section for the three supported
# patterns (network license server / individual license file /
# license-free CodeOcean evaluation).
#
# Build context must be the NNV repo root so `COPY . /home/matlab/nnv`
# includes the engine and tests.
#
# Build:
#   docker build -t nnv3.0 -f code/nnv/examples/NNV3.0/Dockerfile .
#
# Adapted from the official Mathworks reference image:
# https://github.com/mathworks-ref-arch/matlab-dockerfile/blob/main/Dockerfile

# Specify MATLAB release
ARG MATLAB_RELEASE=R2024b

# Specify the list of products to install into MATLAB.
ARG MATLAB_PRODUCT_LIST="MATLAB Computer_Vision_Toolbox Control_System_Toolbox Deep_Learning_Toolbox Image_Processing_Toolbox Optimization_Toolbox Parallel_Computing_Toolbox Statistics_and_Machine_Learning_Toolbox Symbolic_Math_Toolbox System_Identification_Toolbox Deep_Learning_Toolbox_Converter_for_ONNX_Model_Format"

# Specify MATLAB Install Location.
ARG MATLAB_INSTALL_LOCATION="/opt/matlab/${MATLAB_RELEASE}"

# To check the available matlab-deps images, see: https://hub.docker.com/r/mathworks/matlab-deps
FROM mathworks/matlab-deps:${MATLAB_RELEASE}

# Declare build arguments to use at the current build stage.
ARG MATLAB_RELEASE
ARG MATLAB_PRODUCT_LIST
ARG MATLAB_INSTALL_LOCATION

# Install mpm dependencies and Python (with shared library for MATLAB).
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install --no-install-recommends --yes \
    wget \
    ca-certificates \
    python3 \
    python3-pip \
    python3-venv \
    libpython3-dev \
    && apt-get clean \
    && apt-get autoremove \
	&& apt-get install -y git \
    && rm -rf /var/lib/apt/lists/*

# Add "matlab" user and grant sudo permission.
RUN adduser --shell /bin/bash --disabled-password --gecos "" matlab \
    && echo "matlab ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/matlab \
    && chmod 0440 /etc/sudoers.d/matlab

# Set user and work directory.
USER matlab
WORKDIR /home/matlab

# Run mpm to install MATLAB in the target location and delete the mpm installation afterwards.
# If mpm fails to install successfully, then print the logfile in the terminal, otherwise clean up.
# Pass in $HOME variable to install support packages into the user's HOME folder.
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm \
    && chmod +x mpm \
    && sudo HOME=${HOME} ./mpm install \
    --release=${MATLAB_RELEASE} \
    --destination=${MATLAB_INSTALL_LOCATION} \
    --products ${MATLAB_PRODUCT_LIST} \
    || (echo "MPM Installation Failure. See below for more information:" && cat /tmp/mathworks_root.log && false) \
    && sudo rm -rf mpm /tmp/mathworks_root.log \
    && sudo ln -s ${MATLAB_INSTALL_LOCATION}/bin/matlab /usr/local/bin/matlab


# NOTE: MLM_LICENSE_FILE is intentionally NOT baked in here. Supply it
# at `docker run` time via either `-e MLM_LICENSE_FILE=port@server` or
# by mounting a `network.lic` to `/opt/matlab/${MATLAB_RELEASE}/licenses/`.
# See the NNV3.0/README.md "License setup" section.

# Install nnv

COPY --chown=matlab:matlab . /home/matlab/nnv
WORKDIR /home/matlab/nnv

# Create Python virtual environment at .venv (required by ProbVer's cp-star).
RUN python3 -m venv /home/matlab/nnv/.venv \
    && /home/matlab/nnv/.venv/bin/pip install --no-cache-dir -r requirement.txt

ENV PATH="/home/matlab/nnv/.venv/bin:$PATH"
ENV VIRTUAL_ENV="/home/matlab/nnv/.venv"

# NOTE: We deliberately do NOT run `install.m` (which would launch MATLAB
# and require a license at build time). Each NNV3.0 example self-adds NNV
# to the path at runtime via `addpath(genpath(<nnvDir>))`.
