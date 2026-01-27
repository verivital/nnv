# NNV 3.0 CodeOcean Capsule

This capsule runs the NNV 3.0 verification test suite, including:

1. **FairNNV** - Fairness verification on Adult Census dataset
2. **ProbVer** - Probabilistic verification using CP-Star on YOLO benchmark
3. **ModelStar** - Weight perturbation verification on MNIST
4. **VideoStar** - Video classification verification
5. **GNNV** - Graph Neural Network verification on IEEE 24-bus Power Flow

## Model Loading

Different tests use different model loading approaches:

| Test | Model Format | Loading Method |
|------|-------------|----------------|
| **ProbVer** | ONNX (`.onnx`) | `importNetworkFromONNX()` at runtime |
| **FairNNV** | Pre-converted (`.mat`) | Direct `load()` of NNV network |
| **VideoStar** | Pre-converted (`.mat`) | Direct `load()` of NNV network |
| **ModelStar** | Native MATLAB (`.mat`) | `matlab2nnv()` conversion |
| **GNNV** | Native MATLAB (`.mat`) | Direct `load()` of model weights |

### ONNX Runtime Loading (ProbVer)

ProbVer loads ONNX models directly using `importNetworkFromONNX()`. Due to ONNX layer serialization limitations:
- The parallel pool is disabled (`numCores = 1`)
- ONNX-specific layers (`nnet.onnx.layer.*`) are removed from the network
- Falls back to `importONNXNetwork()` if primary loader fails

**ProbVer Python Note:** ProbVer uses CP-Star probabilistic reachability which requires Python (PyTorch, NumPy, SciPy). The Python scripts run via `system()` calls to a virtual environment created during setup. Data is exchanged between MATLAB and Python using `.mat` files (via `scipy.io`).

### Pre-converted Models (FairNNV, VideoStar)

FairNNV and VideoStar use pre-converted `.mat` files containing serialized NNV networks. This approach:
- Avoids runtime ONNX parsing
- Works without the ONNX support package (stub classes provided in `/code/+nnet/+onnx/+layer/`)
- Enables faster startup

## Setup Instructions

### 1. Create CodeOcean Capsule

1. Go to [CodeOcean](https://codeocean.com) and create a new capsule
2. Select **MATLAB R2024b** as the base environment

### 2. Upload Environment Files

Upload the contents of `environment/` to the CodeOcean environment section:

- `Dockerfile` - Custom Docker configuration
- `postInstall` - Post-installation script

### 3. Upload Code Files

Upload all files from `code/` to the CodeOcean `/code` directory:

- `run.m` - Main entry point
- `run_fairnnv.m` - FairNNV test runner
- `run_probver.m` - ProbVer test runner
- `run_modelstar.m` - ModelStar test runner
- `run_videostar.m` - VideoStar test runner
- `run_gnnv.m` - GNNV test runner

### 4. Upload Data Files

Upload the following data files to `/data` on CodeOcean:

```
/data/
├── ICAIF24/                        # FairNNV data
│   ├── adult_onnx/                 # ONNX models
│   │   ├── AC-1.onnx
│   │   └── AC-3.onnx
│   └── data/
│       └── adult_data.mat
├── ProbVer/                        # ProbVer data
│   └── yolo_2023/
│       ├── onnx/                   # ONNX model
│       │   └── TinyYOLO.onnx
│       ├── vnnlib/
│       │   └── *.vnnlib.gz
│       └── instances.csv
├── FORMALISE2025/                  # VideoStar data
│   ├── models/                     # ONNX model
│   │   └── zoomin_4f.onnx
│   └── data/
│       └── ZoomIn/
│           └── *.npy
├── GNNV/                           # GNNV data
│   └── models/
│       ├── gcn_ieee24.mat
│       └── gine_ieee24.mat
└── MNIST/                          # ModelStar data
    └── mnist_model_fc.mat
```

**Source locations in this repository:**

| CodeOcean Path | Source in Repository |
|---------------|---------------------|
| `/data/ICAIF24/adult_onnx/` | `examples/Submission/ICAIF24/onnx/` |
| `/data/ICAIF24/data/` | `examples/Submission/ICAIF24/data/` |
| `/data/ProbVer/yolo_2023/onnx/` | `examples/NNV3.0/ProbVer/yolo_2023/onnx/` |
| `/data/ProbVer/yolo_2023/vnnlib/` | `examples/NNV3.0/ProbVer/yolo_2023/vnnlib/` |
| `/data/FORMALISE2025/models/` | `examples/NNV3.0/VideoStar/models/` |
| `/data/FORMALISE2025/data/` | `examples/NNV3.0/VideoStar/data/` |
| `/data/GNNV/models/` | `codeocean/data/GNNV/models/` |
| `/data/MNIST/mnist_model_fc.mat` | `examples/Tutorial/NN/MNIST/weightPerturb/mnist_model_fc.mat` |

### 5. Configure MATLAB Toolboxes

The following MATLAB toolboxes are required (installed via postInstall):

- Computer Vision Toolbox
- Control System Toolbox
- Deep Learning Toolbox
- Deep Learning Toolbox Converter for ONNX Model Format
- Image Processing Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox
- System Identification Toolbox

### 6. Run the Capsule

Click **Reproducible Run** in CodeOcean.

## Expected Output

Results are saved to `/results/`:

```
/results/
├── test_summary.txt              # Overall test summary
├── FairNNV/
│   ├── fm26_counterfactual_*.csv
│   ├── fm26_individual_*.csv
│   └── fm26_timing_*.csv
├── ProbVer/
│   └── results_summary.csv
├── ModelStar/
│   ├── modelstar_results.csv
│   └── modelstar_summary.txt
├── VideoStar/
│   ├── eps=1_255.csv
│   ├── eps=2_255.csv
│   └── eps=3_255.csv
└── GNNV/
    └── gnnv_results_*.csv
```

## Test Descriptions

### FairNNV
Tests individual and counterfactual fairness on neural networks trained on the Adult Census dataset. Verifies that protected attributes (e.g., sex) don't unfairly influence predictions.

### ProbVer (CP-Star)
Probabilistic verification using the CP-Star reachability method on TinyYOLO object detection benchmark. Tests property satisfaction with coverage and confidence guarantees.

### ModelStar
Weight perturbation verification on MNIST digit classification. Tests robustness to weight perturbations at various magnitudes.

### VideoStar
Video classification verification on ZoomIn dataset. Tests robustness to input perturbations across multiple frames.

### GNNV
Graph Neural Network verification on IEEE 24-bus Power Flow task. Tests GCN and GINE architectures with various epsilon perturbations.

## Notes

- **ProbVer requires GPU**: Select a GPU-enabled machine in CodeOcean for best performance
- **Runtime**: Full test suite may take 30+ minutes depending on configuration
- **Reduced samples**: Some scripts use reduced sample sizes for quicker execution. Modify the config sections to increase coverage.

## Test Coverage

The CodeOcean tests are configured to provide good coverage while keeping runtime reasonable:

| Test | CodeOcean Configuration | Original Configuration | Coverage |
|------|------------------------|------------------------|----------|
| **FairNNV** | 100 samples × 7 epsilon × 2 models | Same | **100%** |
| **ProbVer** | 3 VNNLIB instances | 3 instances (72 available) | **100%** of original |
| **VideoStar** | 10 samples × 3 epsilon | Same | **100%** |
| **GNNV** | 10 scenarios × 3 epsilon × 3 models | Same | **100%** |
| **ModelStar** | 100 images × 3 layers × 4 fracs | Same | **100%** |

### Increasing Coverage

To run more comprehensive tests, modify the config sections in each `run_*.m` file:

- **ProbVer**: Increase `numSamples` to run more VNNLIB instances (up to 72 available)
- **VideoStar**: Expand `config.sampleIndices` beyond `1:10`

## Implementation Differences

The CodeOcean test runners are adapted versions of the original NNV 3.0 example scripts. This section documents the key differences.

**Important:** All tests properly call NNV engine functions for verification. The inline helper functions are either:
1. **Exact copies** from the original scripts (e.g., `perturbationIF`, `L_inf_attack`)
2. **Data preprocessing** that constructs NNV objects (`Star`, `ImageStar`, `VolumeStar`, `GraphStar`)

### NNV Engine Functions Used

| Test | Key NNV Engine Calls |
|------|---------------------|
| FairNNV | `matlab2nnv()`, `net.reach()`, `verify_specification()`, `Star` |
| ProbVer | `Prob_reach()`, `verify_specification()`, `ImageStar`, `load_vnnlib()` |
| VideoStar | `net.verify_robustness()`, `VolumeStar` |
| GNNV | `GNN.reach()`, `GCNLayer`, `GINELayer`, `GraphStar`, `verify_specification()` |
| ModelStar | `WPutils.*`, `matlab2nnv()`, `load_images_MNIST()` |

### General Changes (All Tests)
- **Function format**: CodeOcean scripts are wrapped in functions; originals are scripts
- **Hardcoded paths**: CodeOcean uses `/data/` and `/results/`; originals use relative paths
- **Self-contained**: CodeOcean scripts include all helper functions inline; originals use external dependencies
- **Reduced samples**: Sample counts reduced for faster CodeOcean execution
- **Warning suppression**: CodeOcean scripts suppress warnings for cleaner output

### FairNNV (`run_fairnnv.m`)

| Aspect | Original | CodeOcean |
|--------|----------|-----------|
| Structure | Calls external `adult_verify_fm26.m` and `plot_fm26_results.m` | All logic inline |
| Model loading | Loads ONNX directly via `importNetworkFromONNX` | Uses pre-converted `.mat` files |
| Figure generation | Generates PNG/PDF figures | Skipped (CSV output only) |
| Helper function | Uses external `perturbationIF` | Includes `perturbationIF` inline |
| Lines | ~150 | ~270 |

### ProbVer (`run_probver.m`)

| Aspect | Original | CodeOcean |
|--------|----------|-----------|
| Model loading | `importNetworkFromONNX` | Same, with ONNX layer removal logic |
| Python integration | Uses `pyenv()` for loading `.npz` | Uses `system()` calls only; saves `.mat` |
| Parallel pool | Enabled | Disabled (ONNX layers not serializable) |
| Instance limit | Full VNN-LIB instances | Limited to 3 instances |
| Lines | ~427 | ~450 |

**Key fix (Python)**: The original used MATLAB's `pyenv()` and `py.numpy.load()` to read Python-generated `.npz` files. This fails on CodeOcean because venv Python lacks `libpython.so`. The fix changes Python to save `.mat` files via `scipy.io.savemat()`, and MATLAB loads them natively.

**ONNX workarounds**: ONNX-imported networks contain layers that cannot be serialized to parallel workers:
- Removes `nnet.onnx.layer.*` layers (e.g., `VerifyBatchSizeLayer`, `CustomOutputLayer`) from the layer graph
- Disables parallel pool (`numCores = 1`) because `ONNXParameters` also fail serialization
- Falls back to `importONNXNetwork()` if `importNetworkFromONNX()` fails

### VideoStar (`run_videostar.m`)

| Aspect | Original | CodeOcean |
|--------|----------|-----------|
| Model loading | Loads ONNX via `importNetworkFromONNX` | Uses pre-converted `.mat` files |
| Verification | Calls external `verifyvideo()` from FORMALISE2025 | Inline `L_inf_attack()` and `verify_robustness()` |
| External deps | Requires `FORMALISE2025/src/vvn/` on path | Self-contained |
| npy-matlab | Requires separate npy-matlab installation | Bundled `readNPY` function |
| Lines | ~255 | ~237 |

### GNNV (`run_gnnv.m`)

| Aspect | Original | CodeOcean |
|--------|----------|-----------|
| Scenarios | 10 (indices 1:100:1000) | 10 (indices 1:100:1000) - **matches original** |
| Figure generation | Calls `generate_cav26_dashboard()` and `generate_latex_table()` | Skipped |
| Logging | Optional quiet mode with diary | Always verbose |
| Results storage | Includes `voltage_bounds` and `bound_widths` for figures | Only verification counts |
| Lines | ~558 | ~462 |

### ModelStar (`run_modelstar.m`)

| Aspect | Original | CodeOcean |
|--------|----------|-----------|
| Architecture | Class-based (`EXPT` class with YAML config) | Simple function |
| Configuration | YAML files for experiment definitions | Inline config struct |
| Layers tested | 3 layers (fc_6, fc_5, fc_4) | 3 layers (fc_6, fc_5, fc_4) - **matches original** |
| Perturbation fracs | 4 per layer | 4 per layer - **matches original** |
| Images | 100 | 100 - **matches original** |
| Results | YAML files with plotting support | CSV and TXT summary |
| Model source | Various MNIST architectures | Single `mnist_model_fc.mat` |

### NNV Engine Modifications

The following NNV engine files were modified for CodeOcean compatibility:

**`Direction_trainer.py`** (`engine/nn/Prob_reach/`):
- Changed `np.savez()` to `scipy.io.savemat()` for `.mat` output
- Enables MATLAB to load directions without Python integration

**`ProbReach_ImageStar.m`** (`engine/nn/Prob_reach/`):
- Removed `pyenv()` call and `py.numpy.load()`
- Uses native MATLAB `load()` for `.mat` files

**`Prob_reach.m`** (`engine/utils/`):
- Removed `pyenv` executable path setup (no longer needed)

## Troubleshooting

### "Model file not found"
Ensure all data files are uploaded to the correct paths in `/data/`. Check that ONNX files are in the correct directories.

### "NNV not found"
The postInstall script should clone NNV automatically. If issues persist, manually add NNV to the MATLAB path.

### GPU errors on ProbVer
ProbVer's CP-Star method uses PyTorch which defaults to CPU mode on CodeOcean. For faster execution, select a GPU-enabled compute environment.

### Python/pyenv errors on ProbVer
ProbVer uses Python scripts for CP-Star training. The scripts run via `system()` calls (not MATLAB's Python integration). Data is exchanged using `.mat` files via `scipy.io.savemat()` and MATLAB's `load()`. If you see Python-related errors:
- Ensure the postInstall script created the venv at `/deps/probver_venv`
- Check that torch, numpy, and scipy are installed in the venv
- The run.m script creates a symlink from `nnvroot()/.venv` to the actual venv location

### ONNX layer serialization warnings
When running ProbVer, you may see warnings about ONNX layers not being serializable to parallel workers. This is expected - the parallel pool is disabled automatically to work around this limitation.
