# Semantic Segmentation Network Robustness Verification

Author: Anomynous

## 1. Installation

### 1-a0. OS: Window 10, Linux, (may work on MacOS, some bugs of MPT solvers have been reported)

For the first two steps, more details are available here:

https://github.com/verivital/nnv#installation

### 1-a1. Clone NNV repository using the recursive option as it has submodules:

`git clone --recursive https://github.com/verivital/nnv.git`

### 1-a2. Install NNV (install.m)
In Matlab, navigate to the `code/nnv/` folder. Execute the `install.m` script, which will set up various dependencies (mostly via tbxmanager). This should also set up the path correctly within Matlab so all dependencies are available.

`install`

https://github.com/verivital/nnv/blob/master/code/nnv/install.m

If Matlab is restarted, to work again, either `install.m` or `startup_nnv.m` must be executed again. Alternatively, one can execute `savepath` to update the path after executing install (but in this case, Matlab may need to have been launched with administrative privilege).

`savepath`

https://github.com/verivital/nnv/blob/master/code/nnv/startup_nnv.m

## 2. Run paper reproducibility

Navigate to NeurIPS2020/ (supplementary data provided by authors).

Execute the `reproduce_NeurIPS2020.m` 

