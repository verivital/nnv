# (Artifact) Robustness Verification of Video Classification Neural Networks

This artifact is used to reproduce the results in _Robustness Verification of Video Classification Neural Networks_.

Included in the artifact are the NNV tool, datasets, and scripts used to produce all of the results in the aforementioned paper. The paper introduces a novel abstract set representation for handling layer types common to video classification neural network architectures (3D convolutional, 3D pooling, etc.). The implementation of this abstract set representation is done within the NNV tool, which results in its inclusion in the artifact. Specifically, there exists a subdirectory (`nnv/code/nnv/examples/Submission/FORMALISE2025`) of the artifact which contains the scripts necessary to using the NNV tool for verification and producing the results in the paper.

All results from _Robustness Verification of Video Classification Neural Networks_ were captured using an `Apple M1 Max 10-core CPU@3.20GHz×10` with 64GB of RAM.

# Requirements

The following resources are required to run this artifact:

- MATLAB 2024a with the NNV tool and `npy-matlab` packages installed and added to the path as well as the following toolboxes installed:
  - Computer Vision
  - Control Systems
  - Deep Learning
  - Image Processing
  - Optimization
  - Parallel Computing
  - Statistics and Machine Learning
  - Symbolic Math
  - System Identification
  - [Deep Learning Toolbox Converter for ONNX Model Format](https://www.mathworks.com/matlabcentral/fileexchange/67296-deep-learning-toolbox-converter-for-onnx-model-format)
- A conda environment with Python v3.11. Install rquirements from `requirements.txt`. Make sure to install the source files.
- The datasets which are available for download here: https://doi.org/10.5281/zenodo.14721214

# Installation

This section describes all of the necessary steps for installing tools, dependencies, etc. needed to reproduce the artifact. For the remainder of these instructions when the `FORMALISE2025` directory is referred to we are really referring to the directory at the following path: `nnv/code/nnv/examples/Submission/FORMALISE2025`.

1. Clone the NNV repository (or download artifact for reviewers accessing this way) and ensure all dependencies have been added to the MATLAB path.

   ```
   # Clone NNV
   git clone --recursive https://github.com/verivital/nnv.git
   ```

   Next, we provide instructions for installing both pieces of software as needed for using the artifact. Note that it will be necessary to perform installations with administrator privileges for both pieces of software so that it is possible to `savepath` in MATLAB after completing their respective installation instructions.

2. Install the NNV tool

   To setup the NNV tool, there are a number of toolboxes that we need to first install. Assuming you have already installed MATLAB R2024a, then you can use the `nnv/install_ubuntu.sh` if you are on an Ubuntu machine to speedup this process. Otherwise, you will need to ensure that the toolboxes mentioned under [Requirements](#requirements) are included in your MATLAB installation.

   Additionally, you must install the following support package [Deep Learning Toolbox Converter for ONNX Model Format](https://www.mathworks.com/matlabcentral/fileexchange/67296-deep-learning-toolbox-converter-for-onnx-model-format).

   > Note: Support packages can be installed in MATLAB's HOME tab > Add-Ons > Get Add-ons, search for the support package using the Add-on Explorer and click on the Install button.

   After completing the initial steps, the next step is to run the `nnv/code/nnv/install.m` script in MATLAB. If using a computer with an Apple silicon CPU, then it may be necessary to remove the following line from NNV's `install.m` script:

   ```matlab
   tbxmanager install lcp hysdel cddmex clpmex glpkmex fourier sedumi;
   ```

   Either comment it out or remove it.

   > Note: if you restart Matlab, rerun either `install.m` or `startup_nnv.m`, which will add the necessary dependencies to the path; you alternatively can run savepath after installation to avoid this step after restarting Matlab, but this may require administrative privileges

   All of these instructions are additionally included under heading `# Installation:` on [NNV's README](https://github.com/verivital/nnv/blob/master/README.md). Note, if using a computer with an Apple silicon CPU, then it may be necessary to remove the following line from NNV's `install.m` script:

3. Install the `npy-matlab` package

   For `npy-matlab`, add and save the path to MATLAB with the following commands from the MATLAB interface:

   ```matlab
   >> addpath('/path/to/npy-matlab/npy-matlab')
   >> savepath
   ```

   For reference, the `npy-matlab` package is included in this artifact inside the `FORMALISE2025` directory.

   Instructions are also available on [`npy-matlab`'s README](https://github.com/kwikteam/npy-matlab/blob/master/README.md).

4. Download the `data.tar.gz` file which contains all variations of datasets needed for reproducing results in the desired file structure from the following link: https://doi.org/10.5281/zenodo.14721214.

   After downloading the file, move it to the `FORMALISE2025` directory and unpack it there (`tar -xzf data.tar.gz` if possible, but can also extract files using whatever file explorer is available on machine) so that the file tree now looks like

   ```pseudo
   FORMALISE2025
   ├── data
   │   ├── ZoomIn
   │   ├── ZoomOut
   │   ├── GTSRB
   │   └── STMNIST
   │   └── data.tar.gz
   ...
   ```

   after which the `data.tar.gz` file can be deleted. Please make sure the organization of the folder is the same as shown above.

5. Create a conda environment with `Python v3.11` and install the requirements from `requirements.txt`.

   For producing the results, the Anaconda distribution was used so its installation instructions are provided here; users can also follow general installation instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/). When installing Anaconda in macOS, a `.pkg` file will be provided. Double-click the `.pkg` file and follow the prompts to install Anaconda. For Linux, a `.sh` script will be provided. Executing the command `bash anaconda-latest-Linux-x86_64.sh` and following the prompts will install Anaconda.

   To activate conda in your terminal environment, you will have to restart it. Then, you can test that the installation was successful with

   ```bash
   conda list
   ```

   We can then create a `Python v3.11` environment with command

   ```bash
   conda create --name <env_name> python=3.11
   ```

   Then, activate the environment with and check its Python version

   ```bash
   conda activate <env_name>
   python --version
   ```

   After the environment is activated, install the Python dependencies from `requirements.txt` along with the source files by running the following commands from the `FORMALISE2025` directory in your terminal

   ```bash
   cd /path/to/nnv/code/examples/Submission/FORMALISE2025

   # install requirements
   pip install -r requirements.txt

   # before installing source files, make sure to navigate to the src directory e.g.
   cd src/ # from FORMALISE2025 directory
   pip install -e .
   ```

   Now all Python dependencies have been installed.

6. Modify the `.env` file to add the path to your NNV and npy-matlab directories (the repositories that were cloned earlier). For the `npy-matlab` repository, you'll want to reference the subfolder in the directory also called `npy-matlab`, e.g. `/some/path/to/npy-matlab/npy-matlab`.

7. With all of these steps done, you are now ready to begin reproducing the artifacts.

# Smoke Test Instructions (~1min)

Instructions to quickly test the software needed to reproduce the artifacts of the paper. If no issues arise during the smoke test, you can safely proceed to reproducing all artifacts as described in the below sections.

1. Open a terminal and navigate to the `FORMALISE2025` directory.
2. Make sure the conda environment with the installed dependencies is activated.
3. Run the following command to enable the smoke test script to be executed and then run it:

```bash
chmod +x run_smoketest.sh && ./run_smoketest.sh
```

4. The smoke test will verify a single sample. If the smoke test is successful, then the following message will be output.

```
**********************************************
              Smoke test complete.
**********************************************
```

# Reproducing a Subset of the Results (~1-2 hours)

Assuming the average runtime for the experiments remains as shown in the paper, then it will take approximately 9-10 days to reproduce the results. For that reason, this set of instructions is for reproducing a subset of the results. More specifically, we reproduce the first row of Table 2 from the paper, e.g. the verification results for the 4-frame variation of the Zoom In dataset with the relax verification algorithm. The results will be output to the console after the script completes. Additionally, this script will generate the reachable output range plots used in Figure 7.

1. Open a terminal and navigate to the `FORMALISE2025` directory.
2. Make sure the conda environment with the installed dependencies is activated.
3. Run the following command to begin reproducing a subset of the artifacts:

```bash
chmod +x run_subset_vvn.sh && ./run_subset_vvn.sh
```

# Reproducing a Subset of the Results pt. 2

There is additionally a script to reproduce a single sample from all variations of the datasets. To run this subset, please perform the following:

1. Open a terminal and navigate to the `FORMALISE25` directory.
2. Make sure the conda environment with the installed dependencies is activated.
3. Run the following command to begin reproducing a subset of the artifacts:

```bash
chmod +x run_single_sample_vvn.sh && ./run_single_sample_vvn.sh
```

# Reproducing the Full Results

This set of instructions describes how to reproduce the full results.

1. Open a terminal and navigate to the `FORMALISE2025` directory.
2. Make sure the conda environment with the installed dependencies is activated.
3. Run the following command to begin reproducing the artifacts:

```bash
chmod +x run_vvn.sh && ./run_vvn.sh
```

# Results

| Artifact | Filepath                                                                  | Description                                                                                                                                                     |
| -------- | ------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Table 2  | `results/table2.txt`                                                      | A table of results (PSRV and average runtime) for all experiments.                                                                                              |
| Fig. 7   | `results/reach_stmnist_plot.png` and `results/reach_bad_stmnist_plot.png` | Plots of the reachable output ranges on a sample of STMNIST with different epsilon values (one robust, the other non-robust).                                   |
| Fig. 8   | `results/runtime_comp.png`                                                | A visualization of how average runtime changes for the Zoom In, Zoom Out, and GTSRB datasets as the number of frames and magnitude of the perturbation changes. |

<!-- ### requirements.txt -->
<!-- Numpy could not be upgraded from 1.26.4 to 2.0.0 because of some incompatability with onnxruntime. -->
