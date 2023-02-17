# Verification of Recurrent Neural Networks with Star Reachability
Hoang Dung Tran, Sung Woo Choi, Xiaodong Yang, Tomoya Yamaguchi, Bardh Hoxha and Danil Prokhorov

## HSCC2023

## Source Code:
https://github.com/verivital/nnv


## 1. NNV Installation

### 1.1 System requirements

OS: Windows 10, Linux

RAM: at least 32 GB

Matlab version: 2021b or later

Gurobi optimizer: 9.12 version (or later; 10.00 is tested and works)

Note: Different versions of Gurobi may produce slightly different verification results compared to the paper. The verification results of the paper are produced with Gurobi 9.12 version. Different verification results may be computed if MATLAB optimization is used instead of Gurobi.

Detailed installation of NNV is available at:

https://github.com/verivital/nnv#installation

### 1.2 Clone NNV repository using the recursive option as it has submodules:

`git clone --recursive https://github.com/verivital/nnv.git`

### 1.3 Install NNV (install.m)
In Matlab, navigate to the `code/nnv/` folder. Execute the `install.m` script, which will set up various dependencies (mostly via tbxmanager). This should also set up the path correctly within Matlab, so all dependencies are available.

`install`

https://github.com/verivital/nnv/blob/master/code/nnv/install.m

If Matlab is restarted to work again, either `install.m` or `startup_nnv.m` must be executed again. Alternatively, one can execute `savepath` to update the path after executing install (but in this case, Matlab may need to have been launched with administrative privilege).

`savepath`

https://github.com/verivital/nnv/blob/master/code/nnv/startup_nnv.m

### 1.4 Install Gurobi for MATLAB

1) Dowload Gurobi and extract.

Go to https://www.gurobi.com/downloads/ and download the correct version of Gurobi.

      wget https://packages.gurobi.com/10.0/gurobi10.0.0_linux64.tar.gz

https://www.gurobi.com/documentation/10.0/remoteservices/linux_installation.html recommends installing Gurobi `/opt` for a shared installtion.

      mv gurobi10.0.0_linux64.tar.gz ~/opt/
      tar xvfz gurobi_server10.0.0_linux64.tar.gz

Note: You might have to create the ~/opt/ directory using mkdir ~/opt first.

Move into the directory and extract the content.

      cd ~/opt/
      tar -xzvf gurobi10.0.0_linux64.tar.gz
      rm gurobi10.0.0_linux64.tar.gz


2) Setting up the environment variables.

Open the `~/.bashrc` file.

      vim ~/.bashrc

Add the following lines, replacing {PATH_TO_YOUR_HOME} with the _aboslute_ path to your home directory, and save the file:

      export GUROBI_HOME="{PATH_TO_YOUR_HOME}/opt/gurobi1000/linux64"
      export GRB_LICENSE_FILE="{PATH_TO_YOUR_HOME}/gurobi.lic"
      export PATH="${PATH}:${GUROBI_HOME}/bin"
      export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

Note: If you installed Gurobi or the license file into a different directory, you have to adjust the paths in the first two lines.

After saving, reload .bashrc:

      source ~/.bashrc

3) Acquire your license from https://www.gurobi.com/academia/academic-program-and-licenses/

At `~/opt/gurobi1000/linux64/bin` copy the `grbgetkey` line from the site and enter it into a terminal.

4) Setting up Gurobi for MATLAB

https://www.gurobi.com/documentation/10.0/quickstart_linux/matlab_setting_up_grb_for_.html

To set up Gurobi for MATLAB, you need to launch the Gurobi MATLAB setup script, `gurobi_steup.m` located at `~/opt/gurobi1000/linux64/matlab` in MATLAB.

Add the directory below to MATLAB such that `linprog` function will launch based on `gurobi` optimization instead of MATLAB optimization.

`~/opt/gurobi1000/linux64/examples/matlab`
            

## 2. RnnVerify Installation

### 2.1 Clone RnnVerify

https://github.com/yuvaljacoby/RnnVerify.git

### 2.2 Install RnnVerify

Detailed installation of RnnVerify is available at:

https://github.com/yuvaljacoby/RnnVerify


Install Gurobi:
    [Get Gurobi license](https://www.gurobi.com/downloads/gurobi-optimizer-eula/) (free academic license is available) 
    Install the Gurobi [Python Interface](https://www.gurobi.com/documentation/9.0/quickstart_mac/the_grb_python_interface_f.html)

Compile Marabou:

      mkdir build
	  cd build
	  cmake ..
      cmake --build .

Install python dependencies: 

    pip install -r requirements.txt


## 3. Core results

1) Figure 3. Average verification times of different approaches for the small RNN network. Approaches include exact-star, RNNSVerify, approx-star, RnnVerify.

2) Table 1. Exact verification results for the small RNN.

3) Figure 4. The ranges of all outputs in the first output set in 8 counter output sets for $x_1$ with $T_{max} = 20$ (Table).

4) Table 2. Verification results for the RNN using the approximate and relaxed reachability.

5) Figure 5. Verification performance with different attack bount $\epsilon$ on network $\mathcal{N}_{8,0}$.

6) Figure 5. Verification performance with different attack bount $\epsilon$ on network $\mathcal{N}_{8,0}$.

7) Table 4. Full verification results for $\mathcal{N}_{2,0}$.

8) Table 5. Full verification results for $\mathcal{N}_{2,2}$.

9) Table 6. Full verification results for $\mathcal{N}_{4,0}$.

10) Table 7. Full verification results for $\mathcal{N}_{4,2}$.

11) Table 8. Full verification results for $\mathcal{N}N_{4,4}$.

12) Table 9. Full verification results for $\mathcal{N}N_{8,0}$.

## 4. Reproducing the results in the paper

### 4.1 Reprodcing all NNV results by a single run
-	Navigate to `code/nnv/examples/Submission/HSCC2023/
-	Run `reproduce_HSCC2022.m`

### 4.1.1 Reroducing small NNV results (includes Figure 3, Table 1, Figure 4, Table 2, Figure 5)
-	Navigate to `code/nnv/examples/Submission/HSCC2023/
-	Run `reproduce_HSCC2023_small_results.m`

### 4.1.2 Reroducing full NNV results (includes Figure 3, Table 1, Figure 4, Table 4, Table 5, Table 6, Table 7, Table 8, Table 9, Figure 5)
-	Navigate to `code/nnv/examples/Submission/HSCC2023/
-	Run `reproduce_HSCC2023_full_results.m`

### 4.2 Reproducing NNV results:
#### 4.2.1 Figure 3
-	Navigate to `code/nnv/examples/Submission/HSCC2023/small_RNN`
-	Run `verify_small_RNN.m`

#### 4.2.2 Table 1
-	Navigate to `code/nnv/examples/Submission/HSCC2023/small_RNN`
-	Run `verify_small_RNN.m`

#### 4.2.3 Figure 4
-	Navigate to `code/nnv/examples/Submission/HSCC2023/small_RNN`
-	Run `verify_small_RNN.m`

#### 4.2.4 Table 2
##### For verification results for $\mathcal{N}_{2,0}$:
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_2_0`
-	Run `verify_N2_0.m`

##### For verification results for $\mathcal{N}_{2,2}$:
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_2_2`
-	Run `verify_N2_2.m`

##### For verification results for $\mathcal{N}_{4,0}$:
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_4_0`
-	Run `verify_N4_0.m`

##### For verification results for $\mathcal{N}_{4,2}$:
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_4_2`
-	Run `verify_N4_2.m`

##### For verification results for $\mathcal{N}_{4,4}$:
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_4_4`
-	Run `verify_N4_4.m`

##### For verification results for $\mathcal{N}_{8,0}$:
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_8_0`
-	Run `verify_N4_4.m`

#### 4.2.5 Figure 5
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_8_0`
-	Run `epsilon_vs_robustness_and_time.m`

#### 4.2.6 Table 4
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_2_0`
-	Run `verify_N2_0_full.m`

#### 4.2.7 Table 5
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_2_2`
-	Run `verify_N_2_2_full.m`

#### 4.2.8 Table 6
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_4_0`
-	Run `verify_N_4_0_full.m`

#### 4.2.9 Table 7
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_4_2`
-	Run `verify_N_4_2_full.m`

#### 4.2.10 Table 8
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_4_4`
-	Run `verify_N_4_4_full.m`

#### 4.2.11 Table 9
-	Navigate to `code/nnv/examples/Submission/HSCC2023/N_8_0`
-	Run `verify_N_8_0_full.m`

### 4.3 Reproduing RnnVerify results:
#### 4.3.1 Figure 3
-   Navigate to RnnVerify directory
-   Run `PYTHONPATH=. python3 rnn_experiment/algorithms_compare/experiment.py exp 25`

#### 4.3.2 Table 2, 4, 5, 6, 7, 8, and 9
-   Navigate to RnnVerify directory
-   Run `PYTHONPATH=. python3 rnn_experiment/self_compare/experiment.py exp all 20`

#### 4.3.3 Figure 5
-   Place `code/nnv/examples/Submission/HSCC2023/diff_attack_bounds.py` file to RnnVerify directory
-   Navigate to RnnVerify directory
-   Run `PYTHONPATH=. python3 diff_attack_bounds.py`

