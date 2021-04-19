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

1) Figure 4-a: The robustness, sensitivity and robust IoU of MNIST networks vs. the number of attacked pixels

2) Figure 4-b: The robustness, sensitivity and robust IoU of MNIST networks vs. input size

3) Figure 5: Example of pixel-class reachable sets of MNIST networks

4) Figure 6: The robustness, sensitivity, and robust IoU and reach-times of M2NIST networks vs. the number of attacked pixels

5) Figure 7: Example of pixel-class reachable sets of M2NIST networks

6) Figure 8: Verification time of MNIST networks and reachability time of ReLU layers compared with other layers

7) Table 2: verification reduction of MNIST networks with different relaxation strategies and relaxation factor

8) Figure 9: Conservativeness of different relaxation strategies

10) Figure 10: Compare the performance of the area-based relaxed reachability with DeepZ and DeepPoly


## 3. Run paper reproducibility

### 3.1 Reproducing all results by a single run

Navigate to CAV2021/ (supplementary data provided by authors).

Execute the `reproduce_CAV2021.m` 

### 3.2 Reproducing all results seprarately

#### 3.2.1 Figure 4a, Figure 5a, Figure 5b, Figure 8a

Navigate to CAV2021/MNIST

Execute the `compare_mnist_nets_vs_inputsize.m`

#### 3.2.2 Figure 4b, Figure 8b

Navigate to CAV2021/MNIST

Execute the `compare_mnist_nets_vs_num_attackedpixels.m`

#### 3.2.3 Figure 8c

Navigate to CAV2021/MNIST

Execute the `compare_mnist_ReLU_reachTime_vs_others.m`

#### 3.2.4 Table 2

Navigate to CAV2021/MNIST

Execute the `mnist_net1_verifyTime_vs_relaxFactor.m`

Execute the `mnist_net2_verifyTime_vs_relaxFactor.m`

Execute the `mnist_net3_verifyTime_vs_relaxFactor.m`

#### 3.2.5 Figure 6a, Figure 6b, Figure 7a, Figure 7b

Navigate to CAV2021/M2NIST

Execute the `compare_m2nist_nets_vs_num_attackedpixels.m`

#### 3.2.6 Figure 9

Navigate to CAV2021/RelaxStarCompare

Execute the `mnist01_conservativeness_vs_relaxFactor.m`

#### 3.2.7 Figure 10

Navigate to CAV2021/RelaxStarCompare

Execute the `mnist01_relax_star_area_vs_eran.m`

## 4 Memory Note ***

Since we run the analysis on multiple images in parallel, users may run into memory problem if their computers have less than 32 GB RAM. 





