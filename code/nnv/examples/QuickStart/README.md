# NNV Quick Start Guide

*Part of NNV 3.0 - see [README_NNV3_CONTRIBUTIONS.md](../../../../README_NNV3_CONTRIBUTIONS.md) for new features*

Welcome to NNV! This folder contains simple examples to get you started quickly.

## Prerequisites

Before running these examples, ensure NNV is set up:
```matlab
install  % Run this on first run
```

Every time Matlab is restarted if not saved to path:

```matlab
startup_nnv  % Run this first if not already done or saved to path
```

## Files in This Folder

| File | Description |
|------|-------------|
| `test_installation.m` | Verify NNV is installed correctly |
| `simple_verification.m` | A minimal neural network verification example |

## Step 1: Test Your Installation

```matlab
test_installation
```

This script checks that:
- Core NNV classes are available
- Star set operations work
- Neural network layers function correctly

## Step 2: Your First Verification

```matlab
simple_verification
```

This example demonstrates:
- Creating a simple neural network
- Defining an input set using Star sets
- Computing reachable outputs
- Checking safety properties

## Next Steps

After completing these examples, explore:

1. **MNIST Tutorial** - Verify a real digit classifier
   ```matlab
   run('../Tutorial/NN/MNIST/verify_fc.m')
   ```

2. **ACAS Xu** - Industry benchmark for collision avoidance
   ```matlab
   run('../Tutorial/NN/ACAS_Xu/verify_onnx_vnnlib.m')
   ```

3. **Control Systems** - Neural network controlled systems
   ```matlab
   run('../Tutorial/NNCS/ACC/Verification/verify.m')
   ```

## Troubleshooting

If something doesn't work:
```matlab
check_nnv_setup  % Run diagnostic tool
```

## Getting Help

- Documentation: `help Star`, `help NN`
- Issues: [github.com/verivital/nnv/issues](https://github.com/verivital/nnv/issues)
