# Verification of Deep Convolutional Neural Networks Using ImageStars
Hoang-Dung Tran, Stanley Bak, Weiming Xiang, and Taylor T Johnson

## Source Code:
https://github.com/verivital/nnv

## CAV2020 AE Files:
https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/CAV2020_ImageStar

Note: there is also an AE for the companion tool paper in the repository, so please do not be confused by that, those details can be found here if interested:

Tool paper AE: https://github.com/verivital/run_nnv_comparison/blob/master/README_AE.md

Tool paper AE files: https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/CAV2020

## NNV Manual:
A detailed manual for NNV exists, which illustrates how to create new examples, etc.

https://github.com/verivital/nnv/blob/master/docs/manual.pdf

## CodeOcean Capsule (username: see easychair, password: see easychair):
CodeOcean is a cloud-based execution environment we have configured to reproduce the paper results, discussed in more detail below.

https://codeocean.com/capsule/1314285/

## Demo Video
A recent video demonstration of NNV is available, although note that this does not focus on the ImageStar approach and is more an overview of the whole tool:

[![Alt text](https://img.youtube.com/vi/fLcEwPae5C4/0.jpg)](https://www.youtube.com/watch?v=fLcEwPae5C4)

## REPRODUCING THE RESULTS IN THE PAPER

### 1. The paper contains the following computational elements:

The computational results of the paper consist of three main parts.

1. Comparison of the Zonotope, Polytope and ImageStar methods on three MNIST networks.
2. Comparison of the Polytope and ImageStar methods on VGG16 and VGG19.
3. Comparison of the exact and approximate scheme of the ImageStar methods on VGG16 and VGG19.

Each reproducible figure and table is discussed next, with installation instructions and more below. While one can reproduce any figure and table from a script, they may all be reproduced through this script:

https://github.com/verivital/nnv/blob/master/code/nnv/examples/Submission/CAV2020_ImageStar/reproduce_CAV2020_ImageStar.m

So to reproduce all figures and tables, one may simply:
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/` 
- Run `reproduce_CAV2020_ImageStar`

If run standalone, the resulting figures, tables, and log files will be one-level above the root directory where NNV was cloned, in a `logs` folder.

#### MNIST

Notes: When producing Table 1, 2, and 3, in the case that the reviewer runs into out of memory (OOM) errors, we suggest the reviewer to run the short version of the results by using `compare_star_absdom_short.m` for each table instead of `compare_star_absdom.m`.  This script will produce a small version of the full result. The full version will also take ~10 hours runtime (on a powerful computer like the CodeOcean one described shortly).

1. Figure 8. Comparison of output ranges of the small MNIST classification networks using different approaches.
-	Go to `code/nnv/example/Submission/CAV2020_ImageStar/MNIST_NETS/Small`
-	Run `plot_ranges.m`

2. Table 1. Small MNIST CNN verification results.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/MNIST_NETS/Small`
- Run `compare_star_absdom.m`

3. Table 2. Medium MNIST CNN verification results.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/MNIST_NETS/Medium`
- Run `compare_star_absdom.m`

4. Table 3. Large MNIST CNN verification results.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/MNIST_NETS/Large`
- Run `compare_star_absdom.m`

5. Figure 13 in Appendix. Architectures of the small, medium, and large MNIST CNNs.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/MNIST_NETS/Architecture`
- Run `plot_network_architectures.m`

#### VGG16 and VGG19 Comparison of Polytope and ImageStar Methods

1. Table 4, VGG16 part. VGG verification results.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/VGG16/Compare_Polytope_ImageStar/`
- Run `verify_VGG16.m`

2. Table 4, VGG19 part. VGG verification results.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/VGG19/Compare_Polytope_ImageStar/`
- Run `verify_VGG19.m`

#### VGG16 and VGG19 Comparison of Exact and Approximate ImageStars

1. Table 5, VGG16 part. VGG verification results comparing exact and approximate.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/VGG16/Compare_Exact_vs_Approx`
- Run verify_robustness_delta_e_07.m  and verify_robustness_delta_2e_07.m

2. Table 5, VGG19 part. VGG verification results comparing exact and approximate.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/VGG19/Compare_Exact_vs_Approx`
- Run `verify_robustness_delta_e_07.m` and `verify_robustness_delta_2e_07.m`

3. Figure 9. Exact ranges of VGG19 shows that VGG19 correctly classifies the input image as a bell pepper.
- Go to code/nnv/example/Submission/CAV2020_ImageStar/VGG19/Plot_Figures
- Run `plot_vgg19_exact_range.m`

4. Figure 10. A counterexample shows that VGG19 misclassifies the input image as a strawberry instead of a bell pepper.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/VGG19/Plot_Figures`
- Run `plot_vgg19_counter_example.m`

5. Figure 11. Total reachability time of each type of layers in the VGG19, where the max pooling and ReLU layers dominate the total reachability analysis runtime.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/VGG19/Plot_Figures`
- Run `plot_vgg19_reachTime.m`

6. Figure 12. Number of ImageStars in exact analysis increases with input size.
- Go to `code/nnv/example/Submission/CAV2020_ImageStar/VGG19/Plot_Figures`
- Run `plot_vgg19_inputSize_effect.m`

### 2. System Requirements

#### 2a. CodeOcean

We recommend using the CodeOcean capsule described shortly, which will avoid you having to install anything. The configured CodeOcean capsule has Ubuntu 18.04 and Matlab 2019a, executed inside a Docker containerized environment. The underlying hardware is shared cloud infrastructure, specifically an AWS r5d.4xlarge (16-core, 120GB RAM, see: https://help.codeocean.com/en/articles/1120508-what-machine-will-my-code-run-on).

As CodeOcean is shared infrastructure, exact runtimes are not expected to be reproducible due to varying compute load. We provide scripts that automatically generate results like those presented in the paper, but that may differ in magnitude. Additionally, as some results take significant runtime (order of a day of runtime), we have configured the set up to that produce similar results that execute more quickly.

#### 2b. OS

Window 10 or Linux, but either should work (CodeOcean described below is run on Ubuntu). Mac may work, but has not been tested by the developers due to not having any Mac hardware, and some users who have used Mac have reported dependency problems, so we do not recommend it.

#### 2c. Dependencies

Matlab 2019b or later (may work on earlier versions, but untested)

Matlab toolboxes listed here: https://github.com/verivital/nnv#installation

#### 2d. Hardware

The computer needs > 60 GB RAM.  Note that a computer with less RAM cannot be used to reproduce all results in the paper, which is part of why the CodeOcean set up is recommended as it has 120GB RAM.

### 3. CodeOcean Capsule

3a. NNV is configured to run without installation and without dependence upon users having Matlab licenses through CodeOcean, which has an agreement with the MathWorks to provide publicly available infrastructure for reproducing computational research.

3b. The current draft CodeOcean capsule may be accessed here, a login is required (no public sharing by link is available until the capsule is published and has an updated DOI):

https://codeocean.com/capsule/1314285/

Username: see easychair
Password: see easychair

3c. After opening the capsule through the above URL, one can view code, navigate existing reproduced results, etc.

Prior executions of all results in this paper are available in the single execution `Run 7557000` as well as the three `Runs 7518355, 7502688, and 7438990` together. This process takes about 20 hours total, including the time to build the Docker container, set up the tools, etc., which takes a few minutes for NNV. One can navigate the results from any prior execution, so e.g., one can view the tables and figures generated for this paper at:

`Run 7557000\logs` (all results, using `compare_star_absdom_short.m` for MNIST)

`Run 7518355\logs` (second half of VGG16 results and all VGG19 results; continued from `run 7502688` after correcting path error)

`Run 7502688\logs` (MNIST results, first half of VGG16 results; stopped due to a path error that was corrected)

`Run 7438990\logs` (full version of MNIST examples using e.g. `compare_star_absdom.m` instead of `compare_star_absdom_short.m` )

For example, Figure 8 showing the MNIST comparison of reachable states can be seen at:

`Run 7502688\logs\MNIST\figure8_mnist_small.png`

To re-run all computations and reproduce the results, one selects "Reproducible Run," which will run scripts to reproduce the paper results. One first must clone the capsule under the menu Capsule->Duplicate to set up a variant with write permissions. We next explain what this "Reproducible Run" process does. Note that while we have set up the shared user account with sufficient runtime to re-run the results, we suggest reviewers look at the prior runs and logs just discussed due to the long runtime and limited computational time available on CodeOcean. If there are any issues with this or reviewers experience runtime limitations, please let us know and we will try to have the quotas increased. As a simple test, after duplicating the capsule to make it editable, we would recommend modifying the `run_codeocean.m` script to just run an individual figure/table, such as `code/nnv/example/Submission/CAV2020_ImageStar/MNIST_NETS/Small/plot_ranges.m` instead of the full `reproduce_CAV2020_ImageStar.m`. This can be done by replacing the `reproduce_CAV2020_ImageStar` call (https://github.com/verivital/nnv/blob/master/code/run_codeocean.m#L15) with:

```
cd MNIST_NETS/Small/
plot_ranges
```

3d. The CodeOcean execution process starts by building the Docker container, which first executes this Dockerfile (if the container is already cached, it doesn't rebuild, so runtime can be faster):

https://github.com/verivital/nnv/blob/master/environment/Dockerfile

This subsequently runs a post-installation script that installs further dependencies:

https://github.com/verivital/nnv/blob/master/environment/postInstall

Finally, the main entry point that is executed when selecting "Reproducible Run" is this bash script, which launches Matlab and sets up NNV:

https://github.com/verivital/nnv/blob/master/code/run

Within this `run` shell script, the computational artifacts are reproduced through the `run_codeocean.m` Matlab script:

https://github.com/verivital/nnv/blob/master/code/run_codeocean.m

Within this, the paper results are recreated by running the `reproduce_CAV2020_ImageStar.m` script:

https://github.com/verivital/nnv/blob/master/code/run_codeocean.m#L15

Which executes the `reproduce_CAV2020_ImageStar.m` script:

https://github.com/verivital/nnv/blob/master/code/nnv/examples/Submission/CAV2020_ImageStar/reproduce_CAV2020_ImageStar.m

There are a few other commands listed within the run_codeocean file as further examples, as we have not presented all examples, case studies, etc. that NNV has been evaluated on within this paper due to space, and there are significantly more.

- We will update the published CodeOcean capsule with all reproducible results described in this paper to generate a new DOI so others can publicly reproduce. Due to the way CodeOcean works (by generating DOIs for published capsules), we have not yet published the updated capsule and would do so after acceptance.

- A prior publicly accessible CodeOcean capsule that reproduces an old version of the tests in NNV is available here:

https://doi.org/10.24433/CO.1314285.v1 

- There are some restrictions to using CodeOcean beyond the computational time limits, which is why some of the reproducibility is configured as it is, e.g., it does not support git submodules, which NNV depends upon for integration with HyST, CORA, and NNMT. There are workarounds for this that are deployed in the current CodeOcean setup (e.g., through the shell scripts).

### 4. Manual / Standalone Installation

This setup is if you want to delve more into NNV without relying on CodeOcean, but requires a Matlab installation and license. For the first two steps, more details are available here:

https://github.com/verivital/nnv#installation

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

In Matlab, navigate to the CAV2020_ImageStar submission folder at `code/nnv/examples/Submission/CAV2020_ImageStar`.

One can then execute the `reproduce_CAV2020_ImageStar.m` script discussed earlier, or any of the individual scripts discussed earlier (https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/CAV2020_ImageStar#1-the-paper-contains-the-following-computational-elements):

https://github.com/verivital/nnv/blob/master/code/nnv/examples/Submission/CAV2020_ImageStar/reproduce_CAV2020_ImageStar.m
