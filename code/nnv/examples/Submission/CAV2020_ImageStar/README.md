# Verification of Deep Convolutional Neural Networks Using ImageStars
Hoang-Dung Tran, Stanley Bak, Weiming Xiang, and Taylor T Johnson

## Source Code:
https://github.com/verivital/nnv

## CAV2020 AE Files:
https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/CAV2020_ImageStar

Note: there is also an AE for the companion tool paper in the repository, so please do not be confused by that, those details can be found here if interested:

Tool paper AE: https://github.com/verivital/run_nnv_comparison/blob/master/README_AE.md

Tool paper AE files: https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/CAV2020

## CodeOcean Capsule (username: see easychair, password: see easychair):
https://codeocean.com/capsule/1314285/

## REPRODUCING THE RESULTS IN THE PAPER

### 1. The paper contains the following computational elements:



### 2. System Requirements

#### 2a. CodeOcean

We recommend using the CodeOcean capsule described shortly, which will avoid you having to install anything. The configured CodeOcean capsule has Ubuntu 18.04 and Matlab 2019a, executed inside a Docker containerized environment. The underlying hardware is shared cloud infrastructure, specifically an AWS r5d.4xlarge (16-core, 120GB RAM, see: https://help.codeocean.com/en/articles/1120508-what-machine-will-my-code-run-on).

As CodeOcean is shared infrastructure, exact runtimes are not expected to be reproducible due to varying compute load. We provide scripts that automatically generate results like those presented in the paper, but that may differ in magnitude. Additionally, as some results take significant runtime (order of a day of runtime), we have configured the set up to that produce similar results that execute more quickly.

#### 2b. OS

Window 10 or Linux, but either should work (CodeOcean described below is run on Ubuntu). Mac may work, but has not been tested by the developers due to not having any Mac hardware, and some users who have used Mac have reported dependency problems, so we do not recommend it.

#### 2c. Dependencies

Matlab 2018b or later (may work on earlier versions, but untested)

### 3. CodeOcean Capsule

3a. NNV is configured to run without installation and without dependence upon users having Matlab licenses through CodeOcean, which has an agreement with the MathWorks to provide publicly available infrastructure for reproducing computational research.

3b. The current draft CodeOcean capsule may be accessed here, a login is required (no public sharing by link is available until the capsule is published and has an updated DOI):

https://codeocean.com/capsule/1314285/

Username: see easychair
Password: see easychair

3c. After opening the capsule through the above URL, one can view code, navigate existing reproduced results, etc.

Prior executions of all results in this paper are available in Runs 7518355, 7502688, and 7438990. This process takes about 20 hours, including the time to build the Docker container, set up the tools, etc., which takes a few minutes for NNV. One can navigate the results from any prior execution, so e.g., one can view the tables and figures generated for this paper at:

Run 7518355\logs (second half of VGG16 results and all VGG19 results)

Run 7502688\logs (MNIST results, first half of VGG16 results)

Run 7438990\logs (full version of MNIST examples using e.g. compare_star_absdom.m instead of compare_star_absdom_short.m )

For example, Figure 8 showing the MNIST comparison of reachable states can be seen at:

Run 7502688\logs\MNIST\figure8_mnist_small.png

To re-run all computations and reproduce the results, one selects "Reproducible Run," which will run scripts to reproduce the paper results. One first must clone the capsule under the menu Capsule->Duplicate to set up a variant with write permissions. We next explain what this process does.

This starts by building the Docker container, which first executes this Dockerfile (if the container is already cached, it doesn't rebuild, so runtime can be faster):

https://github.com/verivital/nnv/blob/master/environment/Dockerfile

This subsequently runs a post-installation script that installs further dependencies:

https://github.com/verivital/nnv/blob/master/environment/postInstall

Finally, the main entry point that is executed when selecting "Reproducible Run" is this bash script:

https://github.com/verivital/nnv/blob/master/code/run

Within this shell script, the computational artifacts are reproduced through the run_codeocean.m Matlab script:

https://github.com/verivital/nnv/blob/master/code/run_codeocean.m

Within this, the paper results are recreated by running this reproduce_CAV2020_ImageStar.m script:

https://github.com/verivital/nnv/blob/master/code/run_codeocean.m#L15

Which executes this script:

https://github.com/verivital/nnv/blob/master/code/nnv/examples/Submission/CAV2020_ImageStar/reproduce_CAV2020_ImageStar.m

There are a few other commands listed within the run_codeocean file as further examples, as we have not presented all examples, case studies, etc. that NNV has been evaluated on within this paper due to space, and there are significantly more.

- We will update the published CodeOcean capsule with all reproducible results described in this paper to generate a new DOI so others can publicly reproduce. Due to the way CodeOcean works (by generating DOIs for published capsules), we have not yet published the updated capsule and would do so after acceptance.

- A prior publicly accessible CodeOcean capsule that reproduces an old version of the tests in NNV is available here:

https://doi.org/10.24433/CO.1314285.v1 

- There are some restrictions to using CodeOcean beyond the computational time limits, which is why some of the reproducibility is configured as it is, e.g., it does not support git submodules, which NNV depends upon for integration with HyST, CORA, and NNMT. There are workarounds for this that are deployed in the current CodeOcean setup (e.g., through the shell scripts).

### 4. Manual / Standalone Installation

This setup is if you want to delve more into NNV without relying on CodeOcean, but requires a Matlab installation and license.

#### 4b1. Clone NNV repository using the recursive option as it has submodules:

`git clone --recursive https://github.com/verivital/nnv.git`

#### 4b2. Install NNV (install.m)
In Matlab, navigate to the `code/nnv/` folder. Execute the `install.m` script, which will set up various dependencies (mostly via tbxmanager). This should also set up the path correctly within Matlab so all dependencies are available.

`install`

https://github.com/verivital/nnv/blob/master/code/nnv/install.m

If Matlab is restarted, to work again, either `install.m` or `startup_nnv.m` must be executed again. Alternatively, one can execute `savepath` to update the path after executing install (but in this case, Matlab may need to have been launched with administrative privilege).

`savepath`

https://github.com/verivital/nnv/blob/master/code/nnv/startup_nnv.m

#### 4b3. Run paper reproducibility

In Matlab, navigate to the CAV2020_ImageStar submission folder at code/nnv/examples/Submission/CAV2020_ImageStar.

One can then execute the `reproduce_CAV2020_ImageStar.m` script discussed earlier that executes in CodeOcean:

https://github.com/verivital/nnv/blob/master/code/nnv/examples/Submission/CAV2020_ImageStar/reproduce_CAV2020_ImageStar.m
