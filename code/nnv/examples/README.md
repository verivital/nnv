# NNV Examples

This directory contains examples demonstrating NNV's neural network verification capabilities.

## Getting Started

**New to NNV?** Start with the `Tutorial/` folder for step-by-step examples with explanations.

## Directory Guide

| Folder | Description |
|--------|-------------|
| `Tutorial/` | Beginner-friendly examples with detailed comments and exercises |
| `NN/` | Advanced neural network verification examples |
| `NNCS/` | Neural Network Control Systems verification |
| `Transformer/` | Vision transformer verification |
| `Submission/` | Paper reproduction code organized by venue/year |
| `NNV2.0/` | Legacy examples from NNV 2.0 |
| `other/` | Miscellaneous examples |

## Recommended Learning Path

### 1. Basic Neural Network Verification
```matlab
% Start here - verify a simple fully-connected network on MNIST
run('Tutorial/NN/MNIST/verify_fc.m')
```

### 2. Image Classification with CNNs
```matlab
% Verify a CNN on German Traffic Sign Recognition
run('Tutorial/NN/GTSRB/verify_robust_1.m')
```

### 3. ONNX and VNN-LIB Support
```matlab
% Load ONNX models and VNN-LIB specifications
run('Tutorial/NN/ACAS_Xu/verify_onnx_vnnlib.m')
```

### 4. Neural Network Control Systems
```matlab
% Verify a neural network controlled adaptive cruise control
run('Tutorial/NNCS/ACC/verify.m')
```

## Quick Test

Verify your NNV installation is working:
```matlab
% Create a simple Star set
S = Star([0;0], [1;1]);
disp('Star set created successfully!');

% Create and evaluate a simple network
W = [1 0; 0 1];
b = [0; 0];
layer = FullyConnectedLayer(W, b);
output = layer.evaluate([1; 2]);
disp(['Output: [' num2str(output') ']']);
```

## Example Categories

### Tutorial Examples
- `Tutorial/NN/MNIST/` - MNIST digit classification
- `Tutorial/NN/GTSRB/` - Traffic sign recognition
- `Tutorial/NN/ACAS_Xu/` - Aircraft collision avoidance
- `Tutorial/NN/malware/` - Malware detection
- `Tutorial/NNCS/ACC/` - Adaptive cruise control
- `Tutorial/NNCS/Airplane/` - Aircraft control

### Advanced Examples
- `NN/CNN/` - Convolutional neural networks
- `NN/ACAS_Xu/` - ACAS Xu benchmark suite
- `NNCS/ACC/` - Advanced ACC examples
- `NNCS/Airplane/` - Aircraft dynamics
- `Transformer/ViT/` - Vision Transformer verification

### Paper Reproduction
The `Submission/` folder contains code to reproduce results from published papers:
- `ARCH*` - ARCH-COMP competition entries
- `CAV*` - CAV conference papers

## Running Examples

Most examples can be run directly:
```matlab
cd examples/Tutorial/NN/MNIST
verify_fc  % Run the verification script
```

Or from any directory:
```matlab
run('examples/Tutorial/NN/MNIST/verify_fc.m')
```

## Troubleshooting

**Missing toolboxes?** Run `startup_nnv` to see dependency status.

**Need help?** Visit [github.com/verivital/nnv/issues](https://github.com/verivital/nnv/issues)

## See Also

- [NNV Documentation](https://github.com/verivital/nnv)
- [Star Set Reachability](../engine/set/Star.m)
- [Neural Network Class](../engine/nn/NN.m)
