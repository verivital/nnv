# NNV 2.0: The Neural Network Verification Tool
Diego Manzanas Lopez, Sung Woo Choi, Hoang-Dung Tran, and Taylor T. Johnson

### Source Code:
https://github.com/verivital/nnv

### CAV2023 AE Files:
https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/NNV2.0/CAV2023

### Zenodo DOI:
https://doi.org/10.5281/zenodo.7874494

### CodeOcean Capsule:
CodeOcean is a cloud-based execution environment we have configured to reproduce the paper results, discussed in more detail below. 
Published capsule: https://doi.org/10.24433/CO.0803700.v1

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

To facilitate the review process, we also provide a way to reproduce a subset of the results (Neural ODEs, RNNs and Semantic Segmentation results), as well as a smoke test which should only take a few minutes.


### System Requirements

#### a. CodeOcean

We recommend using the CodeOcean capsule described shortly, which will avoid you having to install anything. The configured CodeOcean capsule has Ubuntu 20.04 and Matlab 2022b, executed inside a Docker containerized environment. The underlying hardware is shared cloud infrastructure, specifically an AWS r5d.4xlarge (16-core, 120GB RAM, see: https://help.codeocean.com/en/articles/1120508-what-machine-will-my-code-run-on).

As CodeOcean is shared infrastructure, exact runtimes are not expected to be reproducible due to varying compute load. We provide scripts that automatically generate results like those presented in the paper, but that may differ in magnitude. Additionally, as some results take significant runtime (order of a day of runtime), we have configured the set up to produce a subset of the results that execute more quickly.

#### b. OS

Window 10 or Linux, but either should work (CodeOcean described below is run on Ubuntu). Mac may work, but has not been tested by the developers due to not having any Mac hardware, and some users who have used Mac have reported dependency problems, so we do not recommend it.

#### c. Dependencies

- Matlab 2022b
- Matlab toolboxes and support packages listed here: https://github.com/verivital/nnv#installation

#### d. Hardware

The computer needs > 30 GB RAM.  Note that a computer with less RAM may not be able to reproduce all results in the paper, which is part of why the CodeOcean set up is recommended as it has 120GB RAM.

### Understanding the results

First, we want to acknowledge a mistake (human error) found after we submitted the manuscript. In Table 2, there is a mismatch in the number of UNSAT instances of property 5 using the relax 25% method. While the manuscript says there are __5__ UNSAT instances, the actual number of UNSAT instances is __6__. We apologize for this mistake, and we will update the final version of the manuscript (if the paper gets accepted), but this error does not have a large impact on the paper. The goal of Table 2 is to showcase the different Star methods available in NNV, and show:
 1) For more precise methods we can verify a larger number of instances.
 2) This comes at a cost as more precise methods are also more expensive (time consuming) to compute.


Now we are ready to look into the generated results. When reproducing the results of the paper (https://github.com/verivital/nnv/blob/master/code/nnv/examples/NNV2.0/Submission/CAV2023/run_cav23.m), the following files are generated:
-	`Table2.txt` - Section 4.1, Table 2 (Verification of ACAS Xu properties 3 and 4)
-	`Table3.txt` - Section 4.1, Table 3 (Verification results of the RL, tllverify and oval21 benchmarks)
-	`results_4.2-4.4.txt` - Verification and computation time results of section 4.2, 4.3 and 4.4
-	`rnn_verification_time.pdf` - Fig. 2 (a) (RNN Computation Times)
-	`acc_nl.pdf` - Fig. 2 (b) (Neural ODE, NNCS Nonlinear ACC)
-	`fpa.pdf` - Fig. 2 (c) (Neural ODE FPA)
-	`transposed_results_0.0001_1.pdf` - Fig. 3 (a) (Target Image)
-	`transposed_results_0.0001_3.pdf` - Fig. 3 (b) (Transposed SSNN)
-	`dilated_results_0.0001_3.pdf` - Fig. 3 (c) (Dilated SSNN)

Notes
- Figures are generated in pdf and png format for CodeOcean, only in pdf format when running locally.

In the generated results files, there are two main categories
1) Verification / Reachability results (i.e., _acc_nl.pdf_)
2) Computation times (i.e., _rnn_verification_time.pdf_)

For the first category, the reachability plots as well as the number of SAT and UNSAT properties verified are reproducible, with the one exception  mentioned above.

For the second one (_computation times_), there are a variety of reasons why these numbers may differ from those in the paper:
- Hardware
  - Move vs less powerful CPU cores.
- OS
  - Linux vs MAC vs Windows
- CodeOcean
  - The underlying hardware is shared cloud infrastructure, specifically an AWS r5d.4xlarge (16-core, 120GB RAM, see: https://help.codeocean.com/en/articles/1120508-what-machine-will-my-code-run-on).)
- Parallelization
  - Table 2 (Verification of ACAS Xu properties 3 and 4) uses the exact method with parallel computation (8 cores in the paper). If the hardware where this is executed has a different number of cores used, these results would differ from those in the paper.

For more details about possible timing differences when running MATLAB in different OS and hardware, see: https://www.mathworks.com/support/requirements/choosing-a-computer.html 


## REPRODUCING THE RESULTS IN THE PAPER


### Option 1: CodeOcean Capsule

NNV is configured to run without installation and without dependence upon users having Matlab licenses through CodeOcean, which has an agreement with MathWorks to provide publicly available infrastructure for reproducing computational research.

The published CodeOcean capsule may be accessed here:

https://codeocean.com/capsule/6689683/

After opening the capsule through the above URL, one can view code, navigate existing reproduced results, etc. In order to run it, you need to create an account: https://codeocean.com/signup .

If one signs up with an academic email, the account will be granted 10 hours of computation a month (see: https://help.codeocean.com/en/articles/1053605-quota-measurement-usage-variation-between-runs-and-what-if-i-run-out). 

If one plans to use CodeOcean to reproduce the results, please request extra computation time (total of 20 hours) in the first phase to have enough computation time to run the full set of results (the run in the published capsule took 14 hours and 34 minutes, but it may vary across runs). After signing up, one can contact CodeOcean with this request here: https://codeocean.com/contact/

Prior execution of all results in this paper is available in the single execution of the __Published Result__ ran on _April 26th, 2023 09:57_ under the __logs__ tab. This process took about 15 hours total, including the time to build the Docker container (CodeOcean environment), set up the tools, etc., which takes a few minutes for NNV. One can navigate the results from any prior execution, so e.g., one can view the tables and figures generated for this paper at:

https://codeocean.com/capsule/6689683/

##### Full results

To re-run all computations and reproduce the results, one selects "Reproducible Run," which will run scripts to reproduce the paper results. One first must clone the capsule under the menu Capsule->Duplicate to set up a variant with write permissions. 

We next explain what this "Reproducible Run" process does. Note that while we have set it up with sufficient runtime to re-run a subset of the results (all except Tables 2 and 3), we suggest reviewers look at the prior runs and logs just discussed due to the long runtime and limited computational time available on CodeOcean. If there are any issues with this or reviewers experience runtime limitations, please let us know and we will try to have the quotas increased.

##### Subset results

To reproduce a subset of the results (Sections 4.2 to 4.4), which reduces the computation time from ~14-15 hours to ~2 hours one can replace (uncomment) the `%RE_cav23_short` call in line 18 of `run_codeocean.m` with `RE_cav23_short`, and replace the `run_cav23` call in line 19 of `run_codeocean.m` with `%run_cav23`. Then, select "Reproducible Run". Once it is complete, all the results __except for__ Table2.txt and Table3.txt will be generated. The results can be accessed under the _logs_ tab in the last executed run.

##### Smoke test

As the smoke test, after duplicating the capsule to make it editable, we would recommend modifying the `run_codeocean.m` script to run only two experiments instead of the full `run_cav23`. This can be done by replacing (commenting out) the `run_cav23` call in line 19 of `run_codeocean.m` with `%run_cav23`, and adding the following lines to the script (lines 20-24):

```
cd NNV_vs_MATLAB/rl_benchmarks
verifyAll;
cd ../../NeuralODEs/FPA
fpa_reach;
plot_fpa;
```

Then, select "Reproducible Run". This should only take a few minutes (after the environment is set up). Once it is complete, there is only one expected output under the _logs_ tab in the last executed run:

`fpa.pdf` - Fig. 2 (c) (Neural ODE FPA)

For every CodeOcean run, all results are generated into /results/logs/ directory and can be accessed once the `Reproducible Run` is finished, in the right column of the capsule under the __logs__ tab.

For more information about the results, please see: https://github.com/verivital/nnv/edit/master/code/nnv/examples/NNV2.0/Submission/CAV2023/readme.md#understanding-the-results .


### Option 2: Manual / Standalone Installation

This setup is if you want to delve more into NNV without relying on CodeOcean, but requires a Matlab installation and license. For the first steps (installation and setup), more details are available here:

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
   Parallel Computing
   Statistics and Machine Learning
   Symbolic Math
   System Identification
 ```
- Support packages
```
   Deep Learning Toolbox Converter for ONNX Model Format (https://www.mathworks.com/matlabcentral/fileexchange/67296-deep-learning-toolbox-converter-for-onnx-model-format)
   Deep Learning Toolbox Verification Library (https://www.mathworks.com/matlabcentral/fileexchange/118735-deep-learning-toolbox-verification-library)
```       
   Note: Support packages can be installed in MATLAB's HOME tab > Add-Ons > Get Add-ons, and then search for the support package using the Add-on Explorer and click on the Install button.


####  Install NNV (install.m)
In Matlab, navigate to the `code/nnv/` folder. Execute the `install.m` script, which will set up various dependencies (mostly via tbxmanager). This should also set up the path correctly within Matlab so all dependencies are available.

`install`

https://github.com/verivital/nnv/blob/master/code/nnv/install.m

If Matlab is restarted, to work again, either `install.m` or `startup_nnv.m` must be executed again. Alternatively, one can execute `savepath` to update the path after executing install (but in this case, Matlab may need to have been launched with administrative privilege).

`savepath`

https://github.com/verivital/nnv/blob/master/code/nnv/startup_nnv.m

##### Run paper reproducibility (full and short)

In MATLAB, navigate to the CAV2023 submission folder at `code/nnv/examples/NNV2.0/Submission/CAV2023`.

One can then execute the `run_cav23.m` script for the full paper results (~15 hours) or `RE_cav23_short.m`, which reproduces all results from Sections 4.2 to 4.4 (~2 hours).

All the tables and figures will be generated in the same folder (nnv/code/nnv/examples/NNV2.0/Submission/CAV2023/). 

Notes
- Tables 2 and 3 are not generated when the short script ("RE_cav23_short.m") is executed.

##### Smoke test

For a smoke test, one can navigate to the folder __nnv/code/nnv/examples/NNV2.0/Submission/CAV2023/__, and copy the following lines in MATLAB's command window:
```
cd NNV_vs_MATLAB/rl_benchmarks
verifyAll;
cd ../../NeuralODEs/FPA
fpa_reach;
plot_fpa;
```

Once it is finished (this should only take a few minutes), one can see one new figure in the folder __nnv/code/nnv/examples/NNV2.0/Submission/CAV2023/__:

`fpa.pdf` - Fig. 2 (c) (Neural ODE FPA)

### Options 3: Docker

We make use of the environment used for the CodeOcean by exporting it to one's local machine. We will describe the most important steps below, but for more information about this process, please see: https://help.codeocean.com/en/articles/2199842-exporting-capsules-and-reproducing-results-on-your-local-machine . This option requires a valid MATLAB license.

__Export capsule__

Navigate to the top left corner on the CodeOcean capsule (https://doi.org/10.24433/CO.0803700.v1) and click on `Capsule` > `Export...` (do not include data, not necessary for this AE) > `Export`.

Extract the files of the zip file downloaded (`capsule-6689683.zip`) in a folder of your choosing.

__Retrieving MATLAB license__

One needs to provide a MATLAB License file `license.lic` to run the code:
- Create a MATLAB License file: 
	For the docker container to run MATLAB, one has to create a new license file for the container.
	Log in with your MATLAB account at https://www.mathworks.com/licensecenter/licenses/
	Click on your license, and then navigate to
	1. "Install and Activate"
	1. "Activate to Retrieve License File"
	1. "Activate a Computer"
	(...may differ depending on how your licensing is set up).
- Choose:
	- Release: `R2022b`
	- Operating System: `Linux`
	- Host ID: `12:34:56:78:9a:bc` # this should match your local machine MAC address
	- Computer Login Name: `root`
	- Activation Label: `<any name>`
- Download the file and place it in the same folder where you have extracted the files.

#### Reproducing results

In your terminal, navigate to the folder where you've extracted the capsule and execute the following command
```shell
cd code
chmod +x run
cd ..
```

##### Full results

Then, execute the following command, adjusting parameters (MAC address) as needed:
```shell
docker run --rm \
  --workdir /code \
  --mac-address=12:34:56:78:9a:bc \ # this should match your local machine MAC address
  --volume "$PWD/license.lic":/MATLAB/licenses/network.lic \
  --volume "$PWD/data":/data \
  --volume "$PWD/code":/code \
  --volume "$PWD/results":/results/logs \
  registry.codeocean.com/published/a8573f9d-86b0-4289-897c-427f78b25d88:v1 ./run
```

If the command above fails, remove the comment about the MAC address and execute it all in one line, e.g.
```shell
docker run --rm --workdir /code --mac-address=12:34:56:78:9a:bc --volume "$PWD/license.lic":/MATLAB/licenses/network.lic --volume "$PWD/data":/data --volume "$PWD/code":/code --volume "$PWD/results":/results/logs registry.codeocean.com/published/a8573f9d-86b0-4289-897c-427f78b25d88:v1 ./run
```

The above command ("docker run ...") will be used for the full results, subset results and smoke test, only difference is the content of the `run_codeocean.m` script. It is very important that the MAC address is changed to your local MAC address, and it must match that of the `license.lic` file previously downloaded (MATLAB license).

##### Subset results

To run a subset of the results requires a modification in one file. To modify the file, navigate to "your_folder/capsule-6689683/code/" and open `run_codeocean.m` using your favorite editor/IDE. 

Then, one can replace (uncomment) the `%RE_cav23_short` call in line 18 of `run_codeocean.m` with `RE_cav23_short`, and replace the `run_cav23` call in line 19 of `run_codeocean.m` with `%run_cav23`. 

Finally, run the same docker command as above (https://github.com/verivital/nnv/edit/master/code/nnv/examples/NNV2.0/Submission/CAV2023/readme.md#full-results-1). Once it is complete, all the results __except for__ Table2.txt and Table3.txt will be generated. The results can be accessed in the folder: "your_folder/capsule-6689683/results/". 

##### Smoke test

Similar to the subset results, we need to modify the `run_codeocean.m` file. To modify the file, navigate to "your_folder/capsule-6689683/code/" and open `run_codeocean.m` using your favorite editor/IDE. 

Then, one can replace (comment out) the `run_cav23` call in line 19 of `run_codeocean.m` with `%run_cav23`, and adding the following lines to the script (lines 20-24):

```
cd NNV_vs_MATLAB/rl_benchmarks
verifyAll;
cd ../../NeuralODEs/FPA
fpa_reach;
plot_fpa;
```

Once it is finished (this should only take a few minutes), one can see one new figure in the folder "your_folder/capsule-6689683/code/":

`fpa.pdf` - Fig. 2 (c) (Neural ODE FPA)

&nbsp;
&nbsp;
&nbsp;

###### Notes
_The hash of the published files in Zenodo and CodeOcean do not fully match those provided here in this GitHub repository due to the inclusion of this `readme.md`. The CodeOcean capsule does not contain this `readme.md` as we created it after we had sent it for publication. The Zenodo files include an older version of this `readme.md`, but we needed to publish the code there first in order to include its DOI in the AE instructions of this `readme.md`._

_Computation time estimates provided are based on paper results and published CodeOcean run._
