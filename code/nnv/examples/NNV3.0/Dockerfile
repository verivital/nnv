# NOTE:
# THIS DOCKERFILE IS A COPY OF THE ORIGINAL NNV DOCKERFILE
# ONLY THIS ONE WILL NOT CLONE NNV BUT INSTEAD MOUNT THE CURRENT
# LOCAL VERSION SO THAT WE CAN TEST CHANGES.

# Building from Dockerfile provided by Mathworks
# https://github.com/mathworks-ref-arch/matlab-dockerfile/blob/main/Dockerfile

# Specify MATLAB release
ARG MATLAB_RELEASE=R2024b

# Specify the list of products to install into MATLAB.
ARG MATLAB_PRODUCT_LIST="MATLAB Computer_Vision_Toolbox Control_System_Toolbox Deep_Learning_Toolbox Image_Processing_Toolbox Optimization_Toolbox Parallel_Computing_Toolbox Statistics_and_Machine_Learning_Toolbox Symbolic_Math_Toolbox System_Identification_Toolbox Deep_Learning_Toolbox_Converter_for_ONNX_Model_Format"

# Specify MATLAB Install Location.
ARG MATLAB_INSTALL_LOCATION="/opt/matlab/${MATLAB_RELEASE}"

# Specify license server information using the format: port@hostname
ARG LICENSE_SERVER="27009@licenseserver.it.vanderbilt.edu"

# To check the available matlab-deps images, see: https://hub.docker.com/r/mathworks/matlab-deps
FROM mathworks/matlab-deps:${MATLAB_RELEASE}

# Declare build arguments to use at the current build stage.
ARG MATLAB_RELEASE
ARG MATLAB_PRODUCT_LIST
ARG MATLAB_INSTALL_LOCATION
ARG LICENSE_SERVER

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


ENV MLM_LICENSE_FILE=$LICENSE_SERVER

# Install nnv

COPY --chown=matlab:matlab . /home/matlab/nnv
WORKDIR /home/matlab/nnv

# Create Python virtual environment at .venv (required by cp_env.m)
RUN python3 -m venv /home/matlab/nnv/.venv \
    && /home/matlab/nnv/.venv/bin/pip install --no-cache-dir -r requirement.txt

ENV PATH="/home/matlab/nnv/.venv/bin:$PATH"
ENV VIRTUAL_ENV="/home/matlab/nnv/.venv"

# Configure MATLAB to use the venv Python and install NNV
RUN matlab -nodisplay -r "pyenv('Version', '/home/matlab/nnv/.venv/bin/python'); cd code/nnv/; install; savepath; exit()"
