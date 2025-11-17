# CLAUDE.md - AI Assistant Guide for NNV

## Project Overview

**NNV (Neural Network Verification)** is a MATLAB toolbox for formal verification of deep neural networks and learning-enabled cyber-physical systems (CPS). It implements reachability-based methods for analyzing neural networks and control systems in autonomous systems applications.

### Key Information
- **Primary Language**: MATLAB (1,221+ .m files)
- **Secondary Languages**: Python (probabilistic reachability), C (performance-critical operations), Java (via HyST)
- **Project Type**: Academic research tool with 40+ publications
- **Main Application Areas**: Neural network verification, safety verification, robustness analysis, autonomous CPS
- **License**: Academic/Research (check with maintainers for commercial use)

## Repository Structure

```
nnv/
├── code/
│   └── nnv/
│       ├── engine/              # Core verification engine
│       │   ├── nn/              # Neural network components
│       │   │   ├── NN.m         # Main network class (1,664 lines)
│       │   │   ├── layers/      # 43 layer implementations
│       │   │   ├── funcs/       # 8 activation functions
│       │   │   └── Prob_reach/  # Probabilistic methods (Python)
│       │   ├── set/             # Set representations (9 files, 6,527 lines)
│       │   │   ├── Star.m       # Star set representation (core)
│       │   │   ├── ImageStar.m  # Image-specific star sets
│       │   │   ├── VolumeStar.m # 3D volume representations
│       │   │   └── Zono.m       # Zonotope representation
│       │   ├── nncs/            # Neural Network Control Systems
│       │   ├── utils/           # Utilities and converters
│       │   ├── cora/            # CORA submodule (reachability)
│       │   ├── hyst/            # HyST submodule (hybrid systems)
│       │   └── matconvnet/      # MatConvNet integration
│       ├── examples/            # Examples and benchmarks
│       │   ├── Tutorial/        # Getting started tutorials
│       │   ├── NN/              # Neural network examples
│       │   ├── NNCS/            # Control systems examples
│       │   └── Submission/      # Publication code (CAV, VNNCOMP, etc.)
│       └── tests/               # Test suite (104 test files)
│           ├── nn/              # NN component tests
│           ├── nncs/            # Control system tests
│           ├── set/             # Set representation tests
│           └── io/              # I/O format tests
├── data/                        # Test datasets (MNIST, GTSRB, etc.)
├── docs/                        # Documentation
├── environment/                 # Environment configuration
├── results/                     # Output directory
└── .github/workflows/           # CI/CD configuration
```

## Core Concepts

### 1. Star Set Reachability (Fundamental Abstraction)

The **Star set** is the core mathematical representation used throughout NNV:

```
x = V * b, where:
  V = [c v1 v2 ... vn]  (basis matrix)
  b = [1 a1 a2 ... an]^T (predicate variables)
  subject to: C*a <= d   (linear constraints)
```

**Key Files**:
- `code/nnv/engine/set/Star.m` - Core star set implementation (2,045 lines)
- `code/nnv/engine/set/ImageStar.m` - Image-specific star sets (1,525 lines)
- `code/nnv/engine/set/Zono.m` - Zonotope representation (606 lines)

### 2. Reachability Methods

NNV supports multiple reachability methods with different precision/performance trade-offs:

- **`'exact-star'`** - Sound and complete (exact) reachability
- **`'approx-star'`** - Sound but over-approximate (faster)
- **`'abs-dom'`** - Abstract domain methods (fastest, most approximate)
- **`'approx-zono'`** - Zonotope-based approximation

### 3. Layer-Based Architecture

Every layer inherits from base classes and implements:
- `evaluate(input)` - Forward propagation
- `reach(input_set, method)` - Reachability analysis
- `reach_star_single_input(input_star)` - Star-based reachability

**Supported Layer Types** (43 total):
- Convolutional: Conv1D, Conv2D, Conv3D, TransposedConv1D, TransposedConv2D
- Activation: ReLU, LeakyReLU, Sigmoid, Tanh, HardSigmoid, Sign, SaturatingLinear
- Pooling: AveragePooling2D, MaxPooling2D, MaxUnpooling2D, GlobalAveragePooling
- Normalization: BatchNormalization, LayerNormalization, ElementwiseAffine
- Recurrent: LSTM, Recurrent
- Special: ODEblock (Neural ODEs), PixelClassification, Softmax

### 4. Neural Network Class (NN.m)

The unified `NN` class (refactored in NNV 2.0) supports all network types:

```matlab
% Key properties:
net.Layers        % Array of layer objects
net.Connections   % Table for DAG networks
net.reachMethod   % Verification method
net.reachSet      % Computed reachable sets per layer

% Supports:
% - FeedForward Neural Networks (FFNNs)
% - Convolutional Neural Networks (CNNs)
% - Semantic Segmentation Networks (SEGNETs)
% - Recurrent Neural Networks (RNNs)
% - Binary Neural Networks (BNNs)
% - DAG networks (ResNet, U-Net, etc.)
```

**Location**: `code/nnv/engine/nn/NN.m` (1,664 lines)

## Development Conventions

### MATLAB Coding Standards

1. **Object-Oriented Design**
   - Extensive use of classes and inheritance
   - Methods organized in class files (not separate functions)
   - Properties typically have explicit access modifiers

2. **Naming Conventions**
   - Classes: PascalCase (e.g., `ImageStar`, `ReluLayer`)
   - Methods: camelCase (e.g., `reach_star_single_input`, `stepReach`)
   - Files match class names exactly (e.g., `Star.m` contains `Star` class)

3. **File Organization**
   - One class per file
   - Test files prefixed with `test_` (e.g., `test_Star.m`)
   - Examples organized by publication or domain

4. **Documentation**
   - Methods include brief header comments
   - Complex algorithms reference publication equations
   - Examples include README.md files

### Git Workflow

**Branch Naming**:
- Feature branches: `claude/claude-md-<session-id>`
- Always develop on designated branch (never push to main without PR)

**Commit Practices**:
- Use descriptive commit messages
- Reference issue numbers when applicable
- Common prefixes: "Fix", "Add", "Update", "Refactor"

**Current Branch**: `claude/claude-md-mi2j7tyn657l85f2-01QJDCAQwKEF4cYwrrnNVGR7`

### Testing Standards

**Location**: `code/nnv/tests/`

**Running Tests**:
```matlab
% Run all tests from tests/ folder:
runtests(pwd, 'IncludeSubfolders', true);

% Run specific test:
runtests('test_Star.m');
```

**Test Organization**:
- Mirror source structure (e.g., `tests/set/star/` for `engine/set/Star.m`)
- Use MATLAB's testing framework
- Some files marked `todo_test.m` for incomplete tests
- CI runs full test suite on push/PR

**Test File Naming**: `test_<ClassName>.m` or `test_<functionality>.m`

## Installation and Dependencies

### Required MATLAB Toolboxes
- Computer Vision Toolbox
- Control Systems Toolbox
- Deep Learning Toolbox
- Image Processing Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox
- System Identification Toolbox

### Support Packages
- **Required**: Deep Learning Toolbox Converter for ONNX Model Format
- **Optional**: VGG16, VGG19, PyTorch converter, TensorFlow converter

### Git Submodules (Critical!)
```bash
# Always clone recursively:
git clone --recursive https://github.com/verivital/nnv.git

# If already cloned:
git submodule update --init --recursive
```

**Submodules**:
1. **CORA** - COntinuous Reachability Analyzer (for plant dynamics)
2. **HyST** - Hybrid Systems Translation (hybrid automata)
3. **NNMT** - Neural Network Model Transformation

### Installation Process
```bash
# Ubuntu automated:
chmod +x install_ubuntu.sh
./install_ubuntu.sh

# Manual (in MATLAB):
cd /path/to/nnv
install  % Runs install.m
```

**What `install.m` does**:
1. Creates `tbxmanager/` directory
2. Installs MPT toolbox + dependencies (GLPK, CDD, Fourier, SeDuMi, etc.)
3. Patches GLPK solver (line 372 adjustment)
4. Runs `startup_nnv.m` to add paths
5. Imports HyST Java libraries

**Important**: After restarting MATLAB, run `startup_nnv.m` (or `savepath` after installation)

### External Solvers
- **GLPK** (GNU Linear Programming Kit) - Primary LP solver
- **linprog** (MATLAB) - Alternative LP solver
- Installed via MPT toolbox dependencies

## Working with NNV

### Common Tasks

#### 1. Loading Neural Networks

```matlab
% From MATLAB formats:
net = NN.parse(matlab_net);  % SeriesNetwork, DAGNetwork, dlnetwork

% From ONNX:
net = importONNXNetwork('model.onnx', 'OutputLayerType', 'classification');
net = NN.parse(net);

% From MAT file:
load('network.mat', 'net');
```

**Converter Utilities**:
- `code/nnv/engine/utils/onnx2nnv.m` - ONNX → NNV
- `code/nnv/engine/utils/matlab2nnv.m` - MATLAB → NNV

#### 2. Defining Input Sets

```matlab
% Star set for single input:
lb = [0.1; 0.2];  % Lower bounds
ub = [0.3; 0.4];  % Upper bounds
input_set = Star(lb, ub);

% ImageStar for images:
img = imread('image.png');
epsilon = 0.01;  % Perturbation bound
input_set = ImageStar(img, epsilon);

% VNNLIB properties:
[lb, ub, prop] = load_vnnlib('property.vnnlib');
```

#### 3. Reachability Analysis

```matlab
% Set reachability method:
net.reachMethod = 'exact-star';  % or 'approx-star', 'abs-dom'

% Compute reachability:
output_set = net.reach(input_set);

% Access layer-wise results:
layer5_output = net.reachSet{5};
```

#### 4. Robustness Verification

```matlab
% Check safety property:
is_safe = Verifier.check_safety(net, input_set, unsafe_region);

% Certified prediction:
[is_robust, counter_example] = verify_robustness_cp(net, input, epsilon, true_label);
```

#### 5. Exporting Results

```matlab
% Export to VNNLIB:
export2vnnlib(net, input_lb, input_ub, 'output.vnnlib');

% Save results:
save('reachability_results.mat', 'output_set', 'net');
```

### Example Workflows

#### Verify MNIST Robustness
```matlab
% 1. Load network
net = NN.parse(trained_network);

% 2. Load image
img = mnist_images(:,:,1);
true_label = mnist_labels(1);

% 3. Define perturbation
epsilon = 0.05;
input_set = ImageStar(img, epsilon);

% 4. Perform reachability
net.reachMethod = 'exact-star';
output_set = net.reach(input_set);

% 5. Check robustness
[is_robust, ~] = verify_robustness_cp(net, img, epsilon, true_label);
```

#### Analyze Control System Safety
```matlab
% 1. Define plant dynamics
plant = NonlinearODE(3, 1, @dynamics_func, reachability_params);

% 2. Load NN controller
controller = NN.parse(controller_net);

% 3. Create NNCS
nncs = NNCS(controller, plant);

% 4. Define initial set
init_set = Star(lb_init, ub_init);

% 5. Compute reachable set
[R, t] = nncs.reach(init_set, time_steps);

% 6. Check safety
is_safe = check_safety(R, unsafe_region);
```

## Key Files Reference

### Core Engine Files
- `code/nnv/engine/nn/NN.m` - Main neural network class (1,664 lines)
- `code/nnv/engine/set/Star.m` - Star set implementation (2,045 lines)
- `code/nnv/engine/set/ImageStar.m` - Image star sets (1,525 lines)
- `code/nnv/engine/utils/Verifier.m` - Safety property verification
- `code/nnv/engine/utils/onnx2nnv.m` - ONNX import
- `code/nnv/engine/utils/load_vnnlib.m` - VNNLIB property loading

### Configuration Files
- `install.m` - Installation script
- `startup_nnv.m` - Path configuration
- `install_ubuntu.sh` - Ubuntu automated install
- `.github/workflows/ci.yml` - CI/CD configuration

### Documentation
- `README.md` - Main documentation
- `code/nnv/examples/Tutorial/readme.md` - Tutorial guide
- `Docker_Instructions.md` - Docker setup
- `docs/README.md` - Course materials

## Common Pitfalls and Best Practices

### 1. Path Management
**Problem**: MATLAB can't find NNV functions after restart
**Solution**: Run `startup_nnv.m` or run `savepath` after installation

### 2. Submodule Issues
**Problem**: Missing CORA, HyST, or NNMT functionality
**Solution**: Ensure submodules are initialized:
```bash
git submodule update --init --recursive
```

### 3. GLPK Solver Errors
**Problem**: LP solver failures during reachability
**Solution**:
- Verify GLPK installation via `tbxmanager`
- Check that install.m patch (line 372) was applied
- Try alternative: `net.reachOption.lp_solver = 'linprog';`

### 4. Memory Issues
**Problem**: Out of memory for large networks
**Solution**:
- Use approximate methods: `net.reachMethod = 'approx-star';`
- Enable parallel computing: `net.reachOption.numCores = 4;`
- Reduce input set size

### 5. Layer Support
**Problem**: Unsupported layer type error
**Solution**:
- Check `code/nnv/engine/nn/layers/` for supported layers (43 types)
- Some layers may need manual conversion or approximation
- Reference: CAV 2023 paper for latest supported layers

### 6. Testing Conventions
**Important**:
- Always run tests before committing major changes
- Tests are located in `code/nnv/tests/`
- Some tests marked `todo_test.m` are incomplete (don't treat as failures)
- CI runs full test suite on push

### 7. MATLAB Version Compatibility
- **Minimum**: MATLAB 2023a
- **Recommended**: MATLAB 2024b
- Older versions may have compatibility issues with newer layers

## Performance Optimization

### Reachability Method Selection
```matlab
% Fastest → Slowest (Most approximate → Most precise)
'abs-dom'       % Abstract domain (fastest, most over-approximate)
'approx-zono'   % Zonotope approximation
'approx-star'   % Approximate star (good balance)
'exact-star'    % Exact star (slowest, most precise)
```

### Parallel Computing
```matlab
% Enable parallel processing:
net.reachOption.numCores = 4;  % Use 4 cores

% For large batch verification:
parfor i = 1:num_images
    results(i) = verify_image(images(i));
end
```

### Layer-wise Approximation
```matlab
% Use different methods per layer:
net.Layers(1).reachMethod = 'exact-star';    % Critical first layer
net.Layers(2:end).reachMethod = 'approx-star';  % Approximate others
```

## Competition and Benchmarking

NNV actively participates in verification competitions:

### VNNCOMP (VNN Competition)
- **Location**: `code/nnv/examples/Submission/VNN_COMP/`
- **Years**: 2021, 2022, 2023, 2024, 2025
- **Focus**: Neural network verification benchmarks

### ARCH-COMP (Applied Verification for Continuous and Hybrid Systems)
- **Location**: `code/nnv/examples/Submission/ARCH_COMP/`
- **Years**: 2019-2024
- **Focus**: NNCS verification benchmarks

**Benchmarking Utilities**:
- VNNLIB format support for standard benchmarks
- Timeout handling and result logging
- Comparison with other tools (e.g., MATLAB's verification library)

## Publications and Examples

### Finding Publication Code
All publication-specific code is in `code/nnv/examples/Submission/`:
- CAV 2020, 2023 - Core NNV papers
- FORMATS 2022 - Neural ODEs
- HSCC 2023 - RNN verification
- FormaliSE 2023, 2024 - Binary NNs, malware detection
- FMICS 2023 - Time series RNNs

### Reproducibility
- Each submission folder includes README with instructions
- Many have CodeOcean capsules (no installation required)
- Git tags for major publications: https://github.com/verivital/nnv/tags

## Tutorials and Learning Resources

**Best Starting Point**: `code/nnv/examples/Tutorial/`

**Key Tutorials**:
1. **MNIST Robustness** - Basic NN verification
2. **GTSRB Traffic Signs** - Image classification verification
3. **Semantic Segmentation** - M2NIST segmentation networks
4. **ACC System** - Adaptive cruise control NNCS verification
5. **Malware Detection** - BODMAS dataset verification

**Documentation Map**: See https://docs.claude.com (for Claude Code-specific features)

## Docker Environment

**Build**:
```bash
docker build . -t nnv
```

**Run**:
```bash
docker run -it nnv
```

**Requirements**: ~20GB disk space

**Dockerfile location**: `./Dockerfile`

## Python Integration

**Location**: `code/nnv/engine/nn/Prob_reach/`

**Purpose**: Probabilistic reachability analysis using deep learning

**Dependencies**:
```
torch
numpy
scipy
```

**Files**:
- `Direction_trainer.py` - Training direction selection
- `Trainer_ReLU.py` - ReLU network training
- `Trainer_Linear.py` - Linear network training

**Setup**: Python virtual environment created by `install_ubuntu.sh`

## Contributing Guidelines

### Before Making Changes
1. Ensure MATLAB and all toolboxes are installed
2. Run `startup_nnv.m` to set up paths
3. Create feature branch: `claude/claude-md-<session-id>-<description>`
4. Run relevant tests before committing

### Making Changes
1. Follow MATLAB OOP conventions
2. Add tests for new functionality (`tests/` folder)
3. Update documentation if adding new features
4. Run full test suite: `runtests('tests', 'IncludeSubfolders', true)`

### Committing
1. Write clear commit messages
2. Reference related issues
3. Don't commit generated files (see `.gitignore`)
4. Push to feature branch (never directly to main)

### Pull Requests
1. Ensure all tests pass
2. Include description of changes
3. Link to related issues
4. Maintainers: Diego Manzanas Lopez, Hoang-Dung Tran

## Ignored Files (.gitignore)

**Never commit**:
- `code/nnv/tbxmanager/` - Auto-generated toolbox manager
- `code/nnv/UUV/` - Unmanned underwater vehicle (legacy)
- `.venv` - Python virtual environment
- `**/*.m~` - MATLAB backup files
- `data/vgg16_cache.mat`, `data/vgg19_cache.mat` - Large cache files
- `**/shuttle_images/` - Test images

## Contact and Support

**Maintainers**:
- Diego Manzanas Lopez (diego.manzanas.lopez@vanderbilt.edu)
- Hoang-Dung Tran (trhoangdung@gmail.com)

**Issues**: https://github.com/verivital/nnv/issues

**Institution**: Vanderbilt University (VERI-VITAL Lab)

**Funding**: AFOSR, DARPA, NSF

## Quick Reference Commands

```matlab
% Installation
install.m                    % Run after cloning

% Startup
startup_nnv.m               % Add paths (run after MATLAB restart)

% Testing
runtests('tests', 'IncludeSubfolders', true)  % All tests
runtests('tests/nn/')                         % NN tests only

% Load network
net = NN.parse(matlab_network);

% Define input set
input_set = Star(lb, ub);

% Reachability
net.reachMethod = 'exact-star';
output_set = net.reach(input_set);

% Verify robustness
[is_robust, counter] = verify_robustness_cp(net, input, epsilon, label);
```

## Version Information

**Current Version**: NNV 2.0 (as of CAV 2023)
**MATLAB Compatibility**: 2023a or newer (2024b recommended)
**Latest Release**: Check https://github.com/verivital/nnv/releases

## Summary for AI Assistants

When working on NNV:
1. **Always** run `startup_nnv.m` in MATLAB before working
2. **Never** commit to main branch directly
3. **Always** run tests before committing (`runtests`)
4. **Remember** this is academic research code - prioritize correctness over speed
5. **Check** publication code in `examples/Submission/` for examples
6. **Use** Star sets as the primary abstraction
7. **Reference** the 40+ publications for algorithmic details
8. **Test** with Tutorial examples before making broad changes
9. **Verify** submodules are initialized (`git submodule update --init --recursive`)
10. **Document** complex changes and reference equation numbers from papers

---

**Last Updated**: 2025-11-17
**Repository**: https://github.com/verivital/nnv
**Primary Citation**: Lopez et al., "NNV 2.0: The Neural Network Verification Tool", CAV 2023
