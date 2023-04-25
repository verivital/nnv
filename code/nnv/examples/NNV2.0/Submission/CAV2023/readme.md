# NNV 2.0: The Neural Network Verification Tool
Diego Manzanas Lopez, Sung Woo Choi, Hoang-Dung Tran, and Taylor T. Johnson

### Source Code:
https://github.com/verivital/nnv

### CAV2023 AE Files:
https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/NNV2.0/CAV2023

### CodeOcean Capsule:
CodeOcean is a cloud-based execution environment we have configured to reproduce the paper results, discussed in more detail below.

TODO: LINK TO PUBLISHED CODEOCEAN CAPSULE

## Overview

We provide 3 ways to reproduce the results in the paper
- CodeOcean (recommended)
- Standalone Installation
- Docker 
 
The paper contains the following computational elements:
- Comparison NNV vs MATLAB (Section 4.1, Tables 2 and 3)
- Neural ODEs (Section 4.2, Fig. 2 (b) and (c))
- RNNs (Section 4.3, Fig. 2(a))
- Semantic Segmentation (Section 4.4, Fig. 3)

To facilitate the review process, we also provide a way to reproduce a subset of the results, which includes Neural ODEs, RNNs and Semantic Segmentation results.


### System Requirements

#### a. CodeOcean

We recommend using the CodeOcean capsule described shortly, which will avoid you having to install anything. The configured CodeOcean capsule has Ubuntu 20.04 and Matlab 2022b (Update 2), executed inside a Docker containerized environment. The underlying hardware is shared cloud infrastructure, specifically an AWS r5d.4xlarge (16-core, 120GB RAM, see: https://help.codeocean.com/en/articles/1120508-what-machine-will-my-code-run-on).

As CodeOcean is shared infrastructure, exact runtimes are not expected to be reproducible due to varying compute load. We provide scripts that automatically generate results like those presented in the paper, but that may differ in magnitude. Additionally, as some results take significant runtime (order of a day of runtime), we have configured the set up to produce a subset of the results that execute more quickly.

#### b. OS

Window 10 or Linux, but either should work (CodeOcean described below is run on Ubuntu). Mac may work, but has not been tested by the developers due to not having any Mac hardware, and some users who have used Mac have reported dependency problems, so we do not recommend it.

#### c. Dependencies

Matlab 2022b

Matlab toolboxes and support packages listed here: https://github.com/verivital/nnv#installation

#### d. Hardware

The computer needs > 16 GB RAM.  Note that a computer with less RAM may not be able to reproduce all results in the paper, which is part of why the CodeOcean set up is recommended as it has 120GB RAM.


## REPRODUCING THE RESULTS IN THE PAPER


### Option 1: CodeOcean Capsule

NNV is configured to run without installation and without dependence upon users having Matlab licenses through CodeOcean, which has an agreement with the MathWorks to provide publicly available infrastructure for reproducing computational research.

The current draft CodeOcean capsule may be accessed here, a login is required (no public sharing by link is available until the capsule is published and has an updated DOI):

https://codeocean.com/capsule/1314285/

Username: xxxx
Password: xxxx

```
DO THIS OR ASK EVERYONE TO CREATE AN ACCOUNT. ACCOUNT HAS A 10 HOURS AS DEFAULT, SO MAY NEED TO ASK FOR EXTRA COMPUTATION
```

After opening the capsule through the above URL, one can view code, navigate existing reproduced results, etc.

Prior executions of all results in this paper are available in the single execution `Run XXXXXXX`. This process takes about 15 hours total, including the time to build the Docker container, set up the tools, etc., which takes a few minutes for NNV. One can navigate the results from any prior execution, so e.g., one can view the tables and figures generated for this paper at:

`Run XXXXXXX\logs` (all results can be found here))

To re-run all computations and reproduce the results, one selects "Reproducible Run," which will run scripts to reproduce the paper results. One first must clone the capsule under the menu Capsule->Duplicate to set up a variant with write permissions. 
We next explain what this "Reproducible Run" process does. Note that while we have set up the shared user account with sufficient runtime to re-run a subset of the results (all except Tables 2 and 3), we suggest reviewers look at the prior runs and logs 
just discussed due to the long runtime and limited computational time available on CodeOcean. If there are any issues with this or reviewers experience runtime limitations, please let us know and we will try to have the quotas increased. As a simple test,
after duplicating the capsule to make it editable, we would recommend modifying the `run_codeocean.m` script to just run an individual script such as `ADD SCRIPT HERE` instead of the full `run_cav23.m`. This can be done by replacing the `run_cav23.m` call (ADD LINK WITH REFERENCE LINE) with:

```
TODO: ADD A COMMAND HERE TO USE AS SMOKE TEST
```

All results are generated into /results/logs/ directory and can be accessed once the `Reproducible Run` is finished, in the right column of the capsule.
-	`Table2.txt` - Section 4.1, Table 2 (Verification of ACAS Xu properties 3 and 4)
-	`Table3.txt` - Section 4.1, Table 3 (Verification results of the RL, tllverify and oval21 benchmarks)
-	`results_4.2-4.4.txt` - Verification and computation time results of section 4.2, 4.3 and 4.4
-	`rnn_verification_time.png` - Fig. 2 (a) (RNN Computation Times)
-	`acc_nl.png` - Fig. 2 (b) (Neural ODE, NNCS Nonlinear ACC)
-	`fpa.png` - Fig. 2 (c) (Neural ODE FPA)
-	`transposed_results_0.0001_1.png` - Fig. 3 (a) (Target Image)
-	`transposed_results_0.0001_3.png` - Fig. 3 (b) (Transposed SSNN)
-	`dilated_results_0.0001_1.png` - Fig. 3 (c) (Dilated SSNN)

Notes
- Tables 2 and 3 are not generated when the short script is executed.
- Figures are generated in pdf and png format.


### Option 2: Manual / Standalone Installation

This setup is if you want to delve more into NNV without relying on CodeOcean, but requires a Matlab installation and license. For the first three steps, more details are available here:

https://github.com/verivital/nnv#installation

####  Clone NNV repository using the recursive option as it has submodules:

`git clone --recursive https://github.com/verivital/nnv.git`

####  Ensure all toolboxes and support packages are installed
- Toolboxes
```
   Computer Vision
   Control Systems
   Deep Learning
   Image Processing
   Optimization
   Parallel Processing
   Symbolic Math
   System Identification
   Statistics and Machine Learning
 ```
- Support packages
```
   Deep Learning Toolbox Converter for ONNX Model Format (https://www.mathworks.com/matlabcentral/fileexchange/67296-deep-learning-toolbox-converter-for-onnx-model-format)
   Deep Learning Toolbox Verification Library (https://www.mathworks.com/matlabcentral/fileexchange/118735-deep-learning-toolbox-verification-library)
```       
   Note: Support packages can be installed in MATLAB's HOME tab > Add-Ons > Get Add-ons, search for the support package using the Add-on Explorer and click on the Install button.


####  Install NNV (install.m)
In Matlab, navigate to the `code/nnv/` folder. Execute the `install.m` script, which will set up various dependencies (mostly via tbxmanager). This should also set up the path correctly within Matlab so all dependencies are available.

`install`

https://github.com/verivital/nnv/blob/master/code/nnv/install.m

If Matlab is restarted, to work again, either `install.m` or `startup_nnv.m` must be executed again. Alternatively, one can execute `savepath` to update the path after executing install (but in this case, Matlab may need to have been launched with administrative privilege).

`savepath`

https://github.com/verivital/nnv/blob/master/code/nnv/startup_nnv.m

#### Run paper reproducibility

In MATLAB, navigate to the CAV2023 submission folder at `code/nnv/examples/NNV2.0/Submission/CAV2023`.

One can then execute the `run_cav23.m` script for the full paper results or `RE_cav23_short.m`, which reproduces all results from Sections 4.2 to 4.4 (about 1-2 hours of computation).

All the tables and figures will be generated in the same folder (nnv/code/nnv/examples/NNV2.0/Submission/CAV2023/). 
-	`Table2.txt` - Section 4.1, Table 2 (Verification of ACAS Xu properties 3 and 4)
-	`Table3.txt` - Section 4.1, Table 3 (Verification results of the RL, tllverify and oval21 benchmarks)
-	`results_4.2-4.4.txt` - Verification and computation time results of section 4.2, 4.3 and 4.4
-	`rnn_verification_time.pdf` - Fig. 2 (a) (RNN Computation Times)
-	`acc_nl.pdf` - Fig. 2 (b) (Neural ODE, NNCS Nonlinear ACC)
-	`fpa.pdf` - Fig. 2 (c) (Neural ODE FPA)
-	`transposed_results_0.0001_1.pdf` - Fig. 3 (a) (Target Image)
-	`transposed_results_0.0001_3.pdf` - Fig. 3 (b) (Transposed SSNN)
-	`dilated_results_0.0001_1.pdf` - Fig. 3 (c) (Dilated SSNN)

Notes
- Tables 2 and 3 are not generated when the short script is executed.

### Docker

To reproduce the results using Docker (no installation required, but MATLAB's license necessary), follow the next steps: 
1.	Go to directory where zip file has been downloaded and extract the files:
```
cd your_folder	  
unzip cav2023_AE_nnv.zip
``` 
2.	Set MAC-address and MATLAB R2022b License
 -	MAC-Address: Substitute the MAC address (now set to 02:42:ac:11:00:02) in the file "run_re" with the MAC address of your local machine.  
 -	MATLAB License: For the docker container to run MATLAB, one must create a new license file for the container. Log in with your MATLAB account at https://www.mathworks.com/licensecenter/licenses/  . Click on your license, and then navigate to  
     - "Install and Activate" > "Activate to Retrieve License File" (...may differ depending on how your licensing is set up).  
     -	Create a new license file with the following data:   
         - MATLAB version: 2022b
         - Host ID (MAC address): MAC address of your local machine (Docker’s default value: 02:42:ac:11:00:02)
         - User: root  
     - Finally, download the generated file "license.lic" and put it into the current directory "your_folder/license/" 

3.	Run the code (from your_folder, where the files were extracted): 
```
chmod +x run_subset
./run_subset	
```
 
The results will be generated in subdirectory “results” (your_folder/results) 
-	`Table2.txt` - Section 4.1, Table 2 (Verification of ACAS Xu properties 3 and 4)
-	`Table3.txt` - Section 4.1, Table 3 (Verification results of the RL, tllverify and oval21 benchmarks)
-	`results_4.2-4.4.txt` - Verification and computation time results of section 4.2, 4.3 and 4.4
-	`rnn_verification_time.pdf` - Fig. 2 (a) (RNN Computation Times)
-	`acc_nl.pdf` - Fig. 2 (b) (Neural ODE, NNCS Nonlinear ACC)
-	`fpa.pdf` - Fig. 2 (c) (Neural ODE FPA)
-	`transposed_results_0.0001_1.pdf` - Fig. 3 (a) (Target Image)
-	`transposed_results_0.0001_3.pdf` - Fig. 3 (b) (Transposed SSNN)
-	`dilated_results_0.0001_1.pdf` - Fig. 3 (c) (Dilated SSNN)

Notes
- Tables 2 and 3 are not generated when the short script is executed.
- Figures are generated in pdf and png format.

## Understanding the results



