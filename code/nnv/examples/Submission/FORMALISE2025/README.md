# (Artifact) Robustness Verification of Video Classification Neural Networks

This artifact is used to reproduce the results in _Robustness Verification of Video Classification Neural Networks_.

All results from _Robustness Verification of Video Classification Neural Networks_ were captured using an `Apple M1 Max 10-core CPU@3.20GHz×10` with 64GB of RAM.

# Requirements

The following resources are required to run this artifact:

- MATLAB 2024a with NNV and `npy-matlab` installed and added to the path.
- A conda environment with Python v3.11. Install rquirements from `requirements.txt`. Make sure to install the source files.

# Installation

This section describes all of the necessary steps for installing tools, dependencies, etc. needed to reproduce the artifact. For the remainder of these instructions when the `FORMALISE2025` directory is referred to we are really referring to the directory at the following path: `nnv/code/nnv/examples/Submission/FORMALISE2025`.

1. Clone NNV and `npy-matlab` and install them ensuring that both have been added to the MATLAB path.

```
# Clone NNV and npy-matlab
git clone --recursive https://github.com/verivital/nnv.git
git clone https://github.com/kwikteam/npy-matlab.git
```

Follow the instructions under heading `# Installation:` on [NNV's README](https://github.com/verivital/nnv/blob/master/README.md). Note, if using a computer with an Apple silicon CPU, then it may be necessary to remove the following line from NNV's `install.m` script:

```matlab
tbxmanager install lcp hysdel cddmex clpmex glpkmex fourier sedumi;
```

Either comment it out or remove it.

For `npy-matlab`, add and save the path to MATLAB. Again, instructions are available on [`npy-matlab`'s README](https://github.com/kwikteam/npy-matlab/blob/master/README.md).

> NOTE: It will be necessary to perform installations with administrator privileges for both of these tools so that it is possible to `savepath` in MATLAB after completing their respective installation instructions.

2. Download the following dataset files from [here](https://drive.google.com/drive/folders/1sXRtSObHLBTeKVss2IA-NGPKljirgD8D?usp=drive_link) into the `FORMALISE2025/data` folder so that the file tree now looks like

```pseudo
FORMALISE2025
├── data
│   ├── ZoomIn
│   ├── ZoomOut
│   ├── GTSRB
│   └── STMNIST
...
```

Note that the `GTSRB` data folder is quite large, so it will likely be downloaded in separate parts. Please make sure the organization of the folder is the same as in the download link and that all filenames match.

3. Create a conda environment with `Python v3.11` and install the requirements from `requirements.txt`. Instructions for installing conda can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) with links under `# Regular installation` to instructions for various operating systems. After confirming that conda has installed successfully by creating a conda environment and making it active, proceed to install the requirements as detailed below. In addition to installing from `requirements.txt`, you must install the source files to the environment. Both of these can be done by running the following commands from the root directory (`FORMALISE2025`):

```bash
# installing requirements
pip install -r requirements.txt

# before installing source files, make sure to navigate to this src directory, e.g.
cd /path/to/FORMALISE2025/src
pip install -e .
```

4. Modify the `.env` file to add the path to your NNV and npy-matlab directories (the repositories that were cloned earlier). For the `npy-matlab` repository, you'll want to reference the subfolder in the directory also called `npy-matlab`, e.g. `/some/path/to/npy-matlab/npy-matlab`.

5. With all of these steps done, you are now ready to begin reproducing the artifacts.

# Smoke Test Instructions

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

# Reproducing a Subset of the Results

Assuming the average runtime for the experiments remains as shown in the paper, then it will take approximately 9-10 days to reproduce the results. For that reason, this set of instructions is for reproducing a subset of the results. More specifically, we reproduce the first row of Table 2 from the paper, e.g. the verification results for the 4-frame variation of the Zoom In dataset with the relax verification algorithm. The results will be output to the console after the script completes. Additionally, this script will generate the reachable output range plots used in Figure 7.

1. Open a terminal and navigate to the `FORMALISE2025` directory.
2. Make sure the conda environment with the installed dependencies is activated.
3. Run the following command to begin reproducing a subset of the artifacts:

```bash
chmod +x run_subset_vvn.sh && run_subset_vvn.sh
```

# Reproducing a Subset of the Results pt. 2

There is additionally a script to reproduce a single sample from all variations of the datasets. To run this subset, please perform the following:

1. Open a terminal and navigate to the `FORMALISE25` directory.
2. Make sure the conda environment with the installed dependencies is activated.
3. Run the following command to begin reproducing a subset of the artifacts:

```bash
chmod +x run_single_sample_vvn.sh && run_single_sample_vvn.sh
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
