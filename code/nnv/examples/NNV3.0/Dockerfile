# NNV 3.0 reproducibility image.
#
# Build from the repository root so the entire NNV checkout is included in the
# build context (a `.dockerignore` at the repo root keeps the context lean):
#
#   docker build \
#     -t nnv3.0 \
#     -f code/nnv/examples/NNV3.0/Dockerfile \
#     --build-arg LICENSE_SERVER=<port>@<host> \
#     .
#
# At runtime, MATLAB still needs a working licence. Either set LICENSE_SERVER at
# build time (recommended) or override at run time with
# `-e MLM_LICENSE_FILE=<port>@<host>` or `-v` mount your `network.lic`.
#
# GPU experiments (ProbVer, GNNV, VideoStar) require the host's NVIDIA driver
# and the container runtime's GPU passthrough enabled:
#
#   docker run --gpus all -it nnv3.0

# Specify MATLAB release
ARG MATLAB_RELEASE=R2024b

# Specify the list of products to install into MATLAB.
ARG MATLAB_PRODUCT_LIST="MATLAB Computer_Vision_Toolbox Control_System_Toolbox Deep_Learning_Toolbox Image_Processing_Toolbox Optimization_Toolbox Parallel_Computing_Toolbox Statistics_and_Machine_Learning_Toolbox Symbolic_Math_Toolbox System_Identification_Toolbox Deep_Learning_Toolbox_Converter_for_ONNX_Model_Format"

# Specify MATLAB Install Location.
ARG MATLAB_INSTALL_LOCATION="/opt/matlab/${MATLAB_RELEASE}"

# MATLAB licence source (port@hostname). Empty by default; override with
# --build-arg LICENSE_SERVER=27000@licenses.example.org or set MLM_LICENSE_FILE
# at run time.
ARG LICENSE_SERVER=""

# To check the available matlab-deps images, see: https://hub.docker.com/r/mathworks/matlab-deps
FROM mathworks/matlab-deps:${MATLAB_RELEASE}

# Re-declare so the args are visible in this build stage.
ARG MATLAB_RELEASE
ARG MATLAB_PRODUCT_LIST
ARG MATLAB_INSTALL_LOCATION
ARG LICENSE_SERVER

# Install mpm dependencies, Python, and helpers used by the experiments.
RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install --no-install-recommends --yes \
        wget \
        ca-certificates \
        git \
        python3 \
        python3-pip \
        python3-venv \
        libpython3-dev \
    && apt-get clean \
    && apt-get autoremove \
    && rm -rf /var/lib/apt/lists/*

# Add "matlab" user with passwordless sudo.
RUN adduser --shell /bin/bash --disabled-password --gecos "" matlab \
    && echo "matlab ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/matlab \
    && chmod 0440 /etc/sudoers.d/matlab

USER matlab
WORKDIR /home/matlab

# Install MATLAB via mpm. mpm itself does NOT validate the licence, so the build
# succeeds without a licence server reachable. Licence is consumed at first
# `matlab` invocation.
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm \
    && chmod +x mpm \
    && sudo HOME=${HOME} ./mpm install \
        --release=${MATLAB_RELEASE} \
        --destination=${MATLAB_INSTALL_LOCATION} \
        --products ${MATLAB_PRODUCT_LIST} \
    || (echo "MPM Installation Failure. See below for more information:" && cat /tmp/mathworks_root.log && false) \
    && sudo rm -rf mpm /tmp/mathworks_root.log \
    && sudo ln -s ${MATLAB_INSTALL_LOCATION}/bin/matlab /usr/local/bin/matlab

# If a build-arg was supplied, bake it into the image. Otherwise leave
# MLM_LICENSE_FILE unset and let the user provide it at `docker run` time.
ENV MLM_LICENSE_FILE=${LICENSE_SERVER}

# Copy the (lean, courtesy of /.dockerignore) NNV checkout into the image.
COPY --chown=matlab:matlab . /home/matlab/nnv
WORKDIR /home/matlab/nnv

# Vendor npy-matlab so VideoStar's run_zoomin_4f.m can load .npy data without
# requiring network access at runtime. If the host already has the submodule
# populated (it's tracked at this path as a git submodule and the COPY above
# brought it in), skip the clone — otherwise fetch a shallow copy.
RUN NPY_DST=/home/matlab/nnv/code/nnv/examples/Submission/FORMALISE2025/npy-matlab; \
    if [ -d "$NPY_DST" ] && [ -n "$(ls -A "$NPY_DST" 2>/dev/null)" ]; then \
        echo "[npy-matlab] already populated (likely via submodule), skipping clone"; \
    else \
        rm -rf "$NPY_DST"; \
        git clone --depth=1 https://github.com/kwikteam/npy-matlab "$NPY_DST"; \
    fi

# Python venv for cp_env / Prob_reach (PyTorch + CUDA wheels).
RUN python3 -m venv /home/matlab/nnv/.venv \
    && /home/matlab/nnv/.venv/bin/pip install --no-cache-dir -r requirement.txt
ENV PATH="/home/matlab/nnv/.venv/bin:$PATH"
ENV VIRTUAL_ENV="/home/matlab/nnv/.venv"

# Install NNV (set up MATLAB paths, run install.m). The script attempts to
# fetch MPT3 from www.tbxmanager.com; that mirror is intermittently
# unreachable. Failure of that step is non-fatal — the modern NNV3.0 examples
# do not require MPT3 at runtime — but we surface a clear log line either way.
RUN matlab -nodisplay -batch "\
    pyenv('Version', '/home/matlab/nnv/.venv/bin/python'); \
    cd('/home/matlab/nnv/code/nnv'); \
    try, install; catch ME, fprintf(2, '[install warning] %s\\n', ME.message); end; \
    savepath; \
    fprintf('[install] verifying core NNV...\\n'); \
    cd('/home/matlab/nnv/code/nnv'); check_nnv_setup(); \
    "

# Optional: extract the AI Verification Library (AIVL) Support Package for the
# ToolComparison experiment. The tarball is non-redistributable MathWorks
# code, so it is NOT bundled in this repo. If you have it, drop it at:
#   code/nnv/examples/NNV3.0/ToolComparison/utils/atva26-aivl.tar.gz
# (or set TOOLCOMPARISON_AIVL_TARBALL at build time). The bulk COPY above
# brings it into the image automatically. If the tarball isn't present,
# the script prints a warning and returns; the rest of the suite builds
# normally and only the MW-side rows of ToolComparison will fail at run
# time. See ToolComparison/README.md for the tarball recipe.
RUN matlab -nodisplay -batch "\
    cd('/home/matlab/nnv/code/nnv/examples/NNV3.0/ToolComparison/utils'); \
    try, run('toolbox_install.m'); catch ME, fprintf(2, '[toolbox_install warning] %s\\n', ME.message); end; \
    "

# Where the experiments expect to write their outputs.
RUN mkdir -p /home/matlab/nnv/repeatability_output

# Default to dropping into a shell at the NNV3.0 examples directory.
WORKDIR /home/matlab/nnv/code/nnv/examples/NNV3.0
CMD ["/bin/bash"]
