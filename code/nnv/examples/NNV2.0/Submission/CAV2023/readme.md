# NNV 2.0: The Neural Network Verification Tool
Diego Manzanas Lopez, Sung Woo Choi, Hoang-Dung Tran, and Taylor T. Johnson

## Source Code:
https://github.com/verivital/nnv

## CAV2023 AE Files:
https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/NNV2.0/CAV2023

## CodeOcean Capsule:
CodeOcean is a cloud-based execution environment we have configured to reproduce the paper results, discussed in more detail below.

TODO: LINK TO PUBLISHED CODEOCEAN CAPSULE


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

The computer needs > 30 GB RAM.  Note that a computer with less RAM may not be able to reproduce all results in the paper, which is part of why the CodeOcean set up is recommended as it has 120GB RAM.


## REPRODUCING THE RESULTS IN THE PAPER

###  The paper contains the following computational elements:

#### a. Comparison NNV vs MATLAB

#### b. Neural ODEs

#### c. RNNs

#### d. Semantic Segmentation


### 3. CodeOcean Capsule

3a. NNV is configured to run without installation and without dependence upon users having Matlab licenses through CodeOcean, which has an agreement with the MathWorks to provide publicly available infrastructure for reproducing computational research.

3b. The current draft CodeOcean capsule may be accessed here, a login is required (no public sharing by link is available until the capsule is published and has an updated DOI):

https://codeocean.com/capsule/1314285/

Username: see easychair
Password: see easychair

3c. After opening the capsule through the above URL, one can view code, navigate existing reproduced results, etc.

Prior executions of all results in this paper are available in the single execution `Run 7557000`. This process takes about 15 hours total, including the time to build the Docker container, set up the tools, etc., which takes a few minutes for NNV. One can navigate the results from any prior execution, so e.g., one can view the tables and figures generated for this paper at:

`Run XXXXXXX\logs` (all results can be found here))

To re-run all computations and reproduce the results, one selects "Reproducible Run," which will run scripts to reproduce the paper results. One first must clone the capsule under the menu Capsule->Duplicate to set up a variant with write permissions. 
We next explain what this "Reproducible Run" process does. Note that while we have set up the shared user account with sufficient runtime to re-run a subset of the results (all except Tables 2 and 3), we suggest reviewers look at the prior runs and logs 
just discussed due to the long runtime and limited computational time available on CodeOcean. If there are any issues with this or reviewers experience runtime limitations, please let us know and we will try to have the quotas increased. As a simple test,
after duplicating the capsule to make it editable, we would recommend modifying the `run_codeocean.m` script to just run an individual script such as `ADD SCRIPT HERE` instead of the full `run_cav23.m`. This can be done by replacing 
the `run_cav23.m` call (ADD LINK WITH REFERENCE LINE) with:

```
TODO: ADD A COMMAND HERE TO USE AS SMOKE TEST
```

```
TODO: EITHER UPDATE THE LINKS TO THE DOCKERFILE AND POSTINTALL TO EXPLAIN THEM, OR JUST REMOVE IT
```

3d. The CodeOcean execution process starts by building the Docker container, which first executes this Dockerfile (if the container is already cached, it doesn't rebuild, so runtime can be faster):

https://github.com/verivital/nnv/blob/master/environment/Dockerfile

This subsequently runs a post-installation script that installs further dependencies:

https://github.com/verivital/nnv/blob/master/environment/postInstall

Finally, the main entry point that is executed when selecting "Reproducible Run" is this bash script, which launches Matlab and sets up NNV:

https://github.com/verivital/nnv/blob/master/code/run

Within this `run` shell script, the computational artifacts are reproduced through the `run_codeocean.m` Matlab script:

https://github.com/verivital/nnv/blob/master/code/run_codeocean.m

Within this, the paper results are recreated by running the `run_cav23.m` script:

https://github.com/verivital/nnv/blob/master/code/run_codeocean.m#L15

Which executes the `run_cav23.m` script:

https://github.com/verivital/nnv/blob/master/code/nnv/examples/NNV2.0/Submission/CAV2023/run_cav23.m

- A prior publicly accessible CodeOcean capsule that reproduces an old version of the tests in NNV is available here:

https://doi.org/10.24433/CO.1314285.v1 

- There are some restrictions to using CodeOcean beyond the computational time limits, which is why some of the reproducibility is configured as it is, e.g., it does not support git submodules, which NNV depends upon for integration with HyST, CORA, NNMT, and onnx2nnv. There are workarounds for this that are deployed in the current CodeOcean setup (e.g., through the shell scripts).

### Manual / Standalone Installation

This setup is if you want to delve more into NNV without relying on CodeOcean, but requires a Matlab installation and license. For the first three steps, more details are available here:

https://github.com/verivital/nnv#installation

#### a. Clone NNV repository using the recursive option as it has submodules:

`git clone --recursive https://github.com/verivital/nnv.git`

#### b. Ensure all toolboxes and support packages are installed
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


#### c. Install NNV (install.m)
In Matlab, navigate to the `code/nnv/` folder. Execute the `install.m` script, which will set up various dependencies (mostly via tbxmanager). This should also set up the path correctly within Matlab so all dependencies are available.

`install`

https://github.com/verivital/nnv/blob/master/code/nnv/install.m

If Matlab is restarted, to work again, either `install.m` or `startup_nnv.m` must be executed again. Alternatively, one can execute `savepath` to update the path after executing install (but in this case, Matlab may need to have been launched with administrative privilege).

`savepath`

https://github.com/verivital/nnv/blob/master/code/nnv/startup_nnv.m

#### 4b3. Run paper reproducibility

In Matlab, navigate to the CAV2023 submission folder at `code/nnv/examples/NNV2.0/Submission/CAV2023`.

One can then execute the `run_cav23.m` script for the full paper results or `RE_cav23_short.m`, which reproduces all results from Sections 4.2 to 4.4 (about 1-2 hours of computation).

