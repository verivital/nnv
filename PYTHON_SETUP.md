# Python Setup for Conformal Prediction (CP) Verification

NNV's Conformal Prediction (CP) verification feature requires a Python environment with PyTorch and related dependencies. This guide covers setup for Windows, macOS, and Linux.

## Overview

The CP verification methods use Python for neural network training of surrogate models. MATLAB calls Python scripts located in `code/nnv/engine/nn/Prob_reach/`:
- `Trainer_Linear.py` - Linear surrogate model training
- `Trainer_ReLU.py` - ReLU surrogate model training
- `Direction_trainer.py` - Direction computation

## Quick Setup

### Windows

```cmd
cd <nnv-root-directory>
python -m venv .venv
.venv\Scripts\activate
pip install -r requirement.txt
```

### macOS / Linux

```bash
cd <nnv-root-directory>
python -m venv .venv
source .venv/bin/activate
pip install -r requirement.txt
```

## Verification

After setup, verify in MATLAB:

```matlab
% Check Python environment is detected
check_nnv_setup

% Or test directly
python_path = cp_env()
```

The output should show the path to your Python executable in the `.venv` directory.

## Running CP Verification

Once configured, you can run CP verification:

```matlab
% Example: CP robustness verification
cd(fullfile(nnvroot(), 'code', 'nnv', 'examples', 'NN', 'cifar10'))
CP_verify_robustness
```

## Dependencies

The `requirement.txt` file includes:
- `torch` - PyTorch for neural network training
- `numpy` - Numerical computing
- `scipy` - Scientific computing (for .mat file I/O)

## GPU Support

By default, CP verification uses CPU. To use GPU acceleration:

1. Install CUDA-enabled PyTorch (see [PyTorch installation guide](https://pytorch.org/get-started/locally/))
2. Set the training device in MATLAB:

```matlab
reachOptions.train_device = 'gpu';
result = verify_robustness_cp(net, IS, reachOptions, target, numClasses);
```

## Troubleshooting

### Python not found

If you see "Python virtual environment not found":
1. Ensure you created the venv in the NNV root directory (where `requirement.txt` is located)
2. Verify the `.venv` folder exists
3. On Windows, check that `.venv\Scripts\python.exe` exists
4. On Unix, check that `.venv/bin/python` exists

### Module not found errors

If Python scripts fail with import errors:
```bash
# Activate the virtual environment first
# Windows:
.venv\Scripts\activate
# Unix:
source .venv/bin/activate

# Then reinstall dependencies
pip install -r requirement.txt
```

### MATLAB can't find Python

Ensure the virtual environment is in the NNV root directory (same level as `requirement.txt`). The `cp_env()` function uses `nnvroot()` to locate the Python executable.

## Related Documentation

- [NNV README](README.md) - Main installation guide
- [CP verification example](code/nnv/examples/NN/cifar10/CP_verify_robustness.m)
- [Prob_reach functions](code/nnv/engine/utils/Prob_reach.m)
