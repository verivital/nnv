# Semantic Segmentation Network Robustness Verification

Author: Anomynous

## 1. Installation

### 1.1. System requirements

OS: Window 10, Linux, (may work on MacOS, some bugs of MPT solvers have been reported)

RAM: at least 32 GB 

Matlab version: 2019b or later

For the first two steps, more details are available here:

https://github.com/verivital/nnv#installation

### 1.2. Clone NNV repository using the recursive option as it has submodules:

`git clone --recursive https://github.com/verivital/nnv.git`

### 1.3. Install NNV (install.m)
In Matlab, navigate to the `code/nnv/` folder. Execute the `install.m` script, which will set up various dependencies (mostly via tbxmanager). This should also set up the path correctly within Matlab so all dependencies are available.

`install`

https://github.com/verivital/nnv/blob/master/code/nnv/install.m

If Matlab is restarted, to work again, either `install.m` or `startup_nnv.m` must be executed again. Alternatively, one can execute `savepath` to update the path after executing install (but in this case, Matlab may need to have been launched with administrative privilege).

`savepath`

https://github.com/verivital/nnv/blob/master/code/nnv/startup_nnv.m

## 2. Cores results

1) Figure 2-a: The robustness, sensitivity and verification time (in second) of MNIST networks vs. the number of attacked pixels

2) Figure 2-b: The robustness, sensitivity and verification time (in second) of MNIST networks vs. input size

3) Figure 3: Example of pixel-class reachable sets of MNIST networks

4) Figure 4: The robustness, sensitivity, verification time (in second) and reach-times of M2NIST networks vs. the number of attacked pixels

5) Figure 5: Example of pixel-class reachable sets of M2NIST networks


## 3. Run paper reproducibility

### 3.1 Reproducing all results by a single run

Navigate to NeurIPS2020/ (supplementary data provided by authors).

Execute the `reproduce_NeurIPS2020.m` 

### 3.2 Reproducing all results seprarately

#### 3.2.1 Figure 2-a

Navigate to NeurIPS2020/MNIST

Execute the `compare_mnist_nets_vs_num_attackedpixels.m`

#### 3.2.2 Figure 2-b and Figure 3

Navigate to NeurIPS2020/MNIST

Execute the `compare_mnist_nets_vs_inputsize.m`

#### 3.2.3 Figure 4 and Figure 5

Navigate to NeurIPS2020/M2NIST

Execute the `compare_m2nist_nets_vs_num_attackedpixels.m`


## 4 Memory Note ***

Since we run the analysis on multiple images in parallel, users may run into memory problem if their computers have less than 32 GB RAM. 





