# (Artifact) Robustness Verification of Video Classification Neural Networks (FormaliSE'25) + Verifying Video Classification Neural Networks for Action Recognition (SoSym)

This artifact is used to reproduce the results in _Robustness Verification of Video Classification Neural Networks_ (FormaliSE'25) and its journal extension _Verifying Video Classification Neural Networks for Action Recognition_ (SoSym).

Included in the artifact are the NNV tool, datasets, and scripts used to produce all of the results. The paper introduces a novel abstract set representation (VideoStars) for handling layer types common to video classification neural network architectures (3D convolutional, 3D pooling, etc.). The implementation is done within the NNV tool, which is included in this artifact.

## Requirements

- **Docker** with the NNV Docker image built from the provided Dockerfile
- **MATLAB R2024b** (included in the Docker image)
- **Python 3.11** with dependencies from `requirements.txt`

## Quick Start (Docker)

### 1. Build the NNV Docker Image

From the root of the NNV repository:

```bash
cd /path/to/nnv
docker build -t nnv .
```

### 2. Start the Docker Container

```bash
docker run -it --name SoSym --mount type=bind,source="$(pwd)",target=/nnv -w /nnv nnv bash
```

### 3. Set Up the Conda Environment (Inside Container)

```bash
# Install Miniconda
wget -P /home/matlab/ https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash /home/matlab/Miniconda3-latest-Linux-x86_64.sh
# Restart terminal or source ~/.bashrc

# Add MATLAB to path
echo 'export MATLABROOT="/opt/matlab/R2024b"' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH="$MATLABROOT/bin/glnxa64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"' >> ~/.bashrc
source ~/.bashrc

# Create conda environment
conda create --name SoSym python=3.11
conda activate SoSym

# Install dependencies
pip install -r /nnv/code/nnv/examples/Submission/FORMALISE2025/requirements.txt

# Install source package
cd /nnv/code/nnv/examples/Submission/FORMALISE2025/src
pip install -e .
```

### 4. Add NNV to MATLAB Path

```bash
cd /nnv
matlab -batch "addpath(pwd); addpath(genpath(pwd)); savepath"
```

### 5. Apply Box.m Modification

Copy the modified `Box.m` file to replace the original in NNV:

```bash
cp /nnv/code/nnv/examples/Submission/FORMALISE2025/src/Box.m /nnv/code/nnv/engine/set/Box.m
```

This modification removes the try-catch optimization for generator construction, using only the loop-based approach for better memory efficiency with large video inputs.

## Data Preparation

### Redistributable Datasets (Zoom In, Zoom Out, GTSRB)

These datasets are available for download from: https://doi.org/10.5281/zenodo.14721214

Download `data.tar.gz` and extract to the `FORMALISE2025` directory:

```bash
cd /nnv/code/nnv/examples/Submission/FORMALISE2025
tar -xzf data.tar.gz
```

### Non-Redistributable Datasets

The following datasets require running data preparation notebooks:

#### ST-MNIST

1. Download the original ST-MNIST dataset
2. Run `src/data_prep/STMNIST/STMNISTDataGeneration.ipynb`

Output files should be saved to:
```
data/STMNIST/stmnistvideo_{16,32,64}f_test_data_seq.npy
data/STMNIST/stmnistvideo_{16,32,64}f_test_labels_seq.npy
```

#### UCF11

1. Download the UCF11 dataset (will be downloaded automatically by notebook)
2. Run notebooks in order:
   - `src/data_prep/UCF11/preprocess_UCF11.ipynb`
   - `src/data_prep/UCF11/create_grayscale_UCF11.ipynb`
   - `src/data_prep/UCF11/postprocess_UCF11.ipynb`

Output files should be saved to:
```
data/UCF11/ucf11_grayscale_{8,16,32}f_verification_data.npy
data/UCF11/ucf11_grayscale_{8,16,32}f_verification_labels.npy
```

#### KTH Actions

1. Download the KTH Actions dataset (will be downloaded automatically by notebook)
2. Run notebooks in order:
   - `src/data_prep/KTHActions/preprocess_KTHActions.ipynb`
   - `src/data_prep/KTHActions/create_grayscale_KTHActions.ipynb`
   - `src/data_prep/KTHActions/postprocess_KTHActions.ipynb`

Output files should be saved to:
```
data/KTHActions/kthactions_grayscale_{8,16,32}f_verification_data.npy
data/KTHActions/kthactions_grayscale_{8,16,32}f_verification_labels.npy
```

## Running Experiments

Navigate to the FORMALISE2025 directory and activate the conda environment:

```bash
cd /nnv/code/nnv/examples/Submission/FORMALISE2025
conda activate SoSym
```

### Smoke Test (~1 min)

Quick test to verify the setup works:

```bash
chmod +x run_smoketest.sh && ./run_smoketest.sh
```

### Subset of Results (~1-2 hours)

Run a subset of experiments for validation:

```bash
chmod +x run_subset_vvn.sh && ./run_subset_vvn.sh
```

### Single Sample Test

Run one sample from each dataset/algorithm combination:

```bash
chmod +x run_single_sample_vvn.sh && ./run_single_sample_vvn.sh
```

### Full Results (~9-10 days)

Run all experiments:

```bash
chmod +x run_vvn.sh && ./run_vvn.sh
```

### Running Individual Datasets

You can also run experiments for specific datasets:

```bash
# MNIST Video (Zoom In/Out) - relax and approx algorithms
python src/run.py --algorithm both  # or --algorithm relax, --algorithm approx

# GTSRB - relax and approx algorithms
python src/run_gtsrb.py --algorithm both

# ST-MNIST - relax and approx algorithms
python src/run_stmnist.py --algorithm both

# KTH Actions - relax only
python src/run_kthactions.py

# UCF11 - relax only
python src/run_ucf11.py
```

## Experiment Coverage

### Tables 3 & 4: Relax and Approx Algorithms

| Dataset | Frames | Samples | Algorithms |
|---------|--------|---------|------------|
| Zoom In | 4, 8, 16 | 100 | relax, approx |
| Zoom Out | 4, 8, 16 | 100 | relax, approx |
| GTSRB | 4, 8, 16 | 215 | relax, approx |
| ST-MNIST | 16, 32, 64 | 100 | relax, approx |
| KTH Actions | 8, 16, 32 | 25 | relax only |

### Table 5: UCF11 (Relax Only)

| Frames | Output Channels | Samples |
|--------|-----------------|---------|
| 8 | 8, 16, 32, 64 | 25 |
| 16 | 8, 16, 32 | 25 |
| 32 | 8, 16 | 25 |

## Results

Results are saved to `/tmp/results/` inside the container. The directory structure is:

```
/tmp/results/
├── MNIST/
│   └── <ds_type>/<algorithm>/<frames>/  # e.g., zoom_in/relax/4/
├── GTSRB/
│   └── gtsrb/<algorithm>/<frames>/
├── STMNIST/
│   └── stmnist/<algorithm>/<frames>/
├── KTHActions/
│   └── kthactions/relax/<frames>/
└── UCF11/
    └── ucf11/relax/<frames>/
```

To copy results to your local machine:

```bash
# From inside the container
cp -r /tmp/results /nnv/code/nnv/examples/Submission/FORMALISE2025/results

# Or from outside the container
docker cp SoSym:/tmp/results ./results
```

| Artifact | Description |
|----------|-------------|
| Table 2/3/4 | PSRV and average runtime for all experiments |
| Table 5 | UCF11 results with varying model sizes |
| Fig. 7 | Reachable output range plots for ST-MNIST |
| Fig. 8 | Average runtime comparison across datasets |

## File Structure

```
FORMALISE2025/
├── .env                    # Docker container paths
├── README.md               # This file
├── requirements.txt        # Python dependencies
├── run_vvn.sh              # Run all experiments
├── run_subset_vvn.sh       # Run subset of experiments
├── run_single_sample_vvn.sh # Run single sample per config
├── run_smoketest.sh        # Quick smoke test
├── data/                   # Dataset files (see Data Preparation)
├── models/                 # Pre-trained ONNX models
├── npy-matlab/             # MATLAB numpy interface
└── src/
    ├── Box.m               # Modified Box.m for NNV
    ├── run.py              # MNIST Video experiments
    ├── run_gtsrb.py        # GTSRB experiments
    ├── run_stmnist.py      # ST-MNIST experiments
    ├── run_ucf11.py        # UCF11 experiments
    ├── run_kthactions.py   # KTH Actions experiments
    ├── vvn/                # Verification library
    ├── analysis/           # Result analysis scripts
    └── data_prep/          # Data preparation notebooks
        ├── KTHActions/     # KTH Actions preprocessing
        ├── UCF11/          # UCF11 preprocessing
        └── STMNIST/        # ST-MNIST preprocessing
```

## Hardware Used for Original Results

All results from the FormaliSE'25 paper were captured using:
- Apple M1 Max 10-core CPU @ 3.20GHz
- 64GB RAM

All results from the SoSym extension were captured using:
- Intel Xeon Gold 6238R Processor @ 2.20GHz
- 512GB RAM

## Troubleshooting

### MATLAB Engine Issues

If you encounter MATLAB engine connection errors:

```bash
# Ensure MATLAB paths are set
export MATLABROOT="/opt/matlab/R2024b"
export LD_LIBRARY_PATH="$MATLABROOT/bin/glnxa64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
```

### Memory Issues

For large video inputs, the modified `Box.m` uses a loop-based generator construction that is more memory efficient. Make sure you've applied the Box.m modification as described in the setup.

### Permission Issues

Results are saved to `/tmp/` to avoid permission issues on mounted volumes. Remember to copy results out before destroying the container.
