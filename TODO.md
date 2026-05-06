# NNV Toolbox - TODO and Improvement Roadmap

**Last Updated**: 2025-11-17
**Source**: Comprehensive codebase analysis excluding Submission folder

---

## Table of Contents

1. [Critical Fixes (Immediate Action Required)](#critical-fixes)
2. [High Priority - Near Term (1-4 weeks)](#high-priority-near-term)
3. [Medium Priority - Short Term (1-3 months)](#medium-priority-short-term)
4. [Long Term Improvements (3-12 months)](#long-term-improvements)
5. [Testing and Quality](#testing-and-quality)
6. [Documentation](#documentation)
7. [Performance Optimizations](#performance-optimizations)
8. [Refactoring Opportunities](#refactoring-opportunities)

---

## Critical Fixes (Immediate Action Required)

### 🔴 CRITICAL BUG - Must Fix Immediately

#### 1. GlobalAveragePooling1DLayer Variable Name Error
- **File**: `code/nnv/engine/nn/layers/GlobalAveragePooling1DLayer.m`
- **Line**: 112
- **Issue**: Uses undefined variable `in_image` instead of `input` in parfor loop
- **Impact**: Runtime error when using this layer
- **Fix**: Change `in_image.V` to `input.V` (and similar references)
- **Estimated Effort**: 5 minutes

```matlab
% Current (BROKEN):
parfor i=1:length(in_image)
    image(i) = Star(in_image.V(index, :), in_image.C, in_image.d);
end

% Should be:
parfor i=1:length(input)
    image(i) = Star(input.V(index, :), input.C, input.d);
end
```

---

## High Priority - Near Term (1-4 weeks)

### User Experience Improvements

#### 2. Add NN.parse() Static Method
- **Priority**: HIGH
- **Effort**: 1 day
- **Issue**: Documentation mentions `NN.parse()` but it doesn't exist
- **Current Workaround**: Users must use `matlab2nnv()` directly
- **Implementation**: Add static method to NN class that wraps matlab2nnv and onnx2nnv

```matlab
classdef NN
    methods(Static)
        function net = parse(network, varargin)
            if ischar(network) || isstring(network)
                % ONNX file path
                net = onnx2nnv(network, varargin{:});
            else
                % MATLAB network object
                net = matlab2nnv(network);
            end
        end
    end
end
```

#### 3. Create ReachOptions Class with Validation
- **Priority**: HIGH
- **Effort**: 2-3 days
- **Issue**: Unstructured options lead to errors and no autocomplete support
- **Implementation**: Create validated options class

```matlab
classdef ReachOptions
    properties
        reachMethod (1,1) string {mustBeMember(reachMethod, ...
            ["exact-star", "approx-star", "abs-dom", "approx-zono"])} = "approx-star"
        numCores (1,1) double {mustBePositive} = 1
        relaxFactor (1,1) double {mustBeInRange(relaxFactor,0,1)} = 0
        displayProgress (1,1) logical = false
        lpSolver (1,1) string {mustBeMember(lpSolver, ["linprog", "glpk"])} = "linprog"
    end
end
```

#### 4. Improve Error Messages
- **Priority**: HIGH
- **Effort**: 1 week
- **Files**: Throughout codebase (328 error() calls found)
- **Implementation**: Create centralized error helper with suggestions

**Examples to fix**:
- `NN.m:107`: "Invalid number of inputs" → Add "Got X arguments. See 'help NN'"
- `matlab2nnv.m:222`: "Unsupported Class of Layer" → List supported layers
- `load_vnnlib.m:113`: Generic error → Explain VNNLIB format requirements

#### 5. Add Installation Verification Script
- **Priority**: HIGH
- **Effort**: 2 days
- **File**: Create `install_check.m`
- **Purpose**: Verify MATLAB version, toolboxes, GLPK, submodules, Python env

```matlab
function status = install_check()
    % Check MATLAB version (>= 2023a)
    % Check required toolboxes
    % Test GLPK solver
    % Verify submodule initialization
    % Test Python environment (if used)
    % Report results with actionable suggestions
end
```

### Incomplete Implementations

#### 6. Complete UpsampleLayer.reach_single_input
- **Priority**: HIGH
- **Effort**: 3-5 days
- **File**: `code/nnv/engine/nn/layers/UpsampleLayer.m:107`
- **Issue**: Core reachability function not implemented
- **TODO Comment**: "implement this function, just need to modify dimensions of ImageStar or convert Star to ImageStar. Should also support ImageZono and Zono"

#### 7. Complete GlobalAveragePooling Reach Methods
- **Priority**: MEDIUM-HIGH
- **Effort**: 3-5 days each
- **Files**:
  - `GlobalAveragePooling2DLayer.m:101`
  - `GlobalAveragePooling1DLayer.m:99`
- **Issue**: Both have incomplete reach() methods marked TODO

#### 8. Add Tensor Shape Support to ElementwiseAffineLayer
- **Priority**: MEDIUM
- **Effort**: 1 week
- **File**: `code/nnv/engine/nn/layers/ElementwiseAffineLayer.m`
- **Lines**: 98, 125
- **Issue**: Throws errors for unsupported tensor shapes
- **TODO**: "add support for other tensor shapes (Scale/Offset)"

### Example Updates - Critical for User Experience

#### 9. Update All NNCS/ACC Examples (CRITICAL)
- **Priority**: HIGH
- **Effort**: 2-3 days
- **Files**: 10+ files in `code/nnv/examples/NNCS/ACC/Verification/`
- **Issue**: ALL use deprecated `LayerS` and `FFNNS` classes
- **Impact**: High-visibility tutorial example is broken/outdated

**Files to update**:
- Scenarios 1: verify_controller_*.m, falsify_controller_*.m (5 files)
- Scenarios 2: reach_nncACCsystem.m, reach_nncACCsystem_video.m

**Migration Pattern**:
```matlab
# OLD (current):
L = LayerS(weights{1, i}, bias{i, 1}, 'poslin');
Layers = [Layers L];
Controller = FFNNS(Layers);

# NEW (target):
Controller = matlab2nnv(matlab_network);
# OR load pre-trained and use NN.parse()
```

#### 10. Update NN/MNIST Examples
- **Priority**: HIGH
- **Effort**: 1 week
- **Files**: 15+ files in `code/nnv/examples/NN/MNIST/`
- **Issue**: Core demonstration uses old API
- **Impact**: Creates confusion for new users

**Affected directories**:
- Main directory: verify_L_infinity.m and others
- mnist_1_140/, mnist_2_250/, mnist_3_1000/ (12+ files)

#### 11. Update NNCS Examples
- **Priority**: MEDIUM-HIGH
- **Effort**: 1-2 weeks
- **Directories**:
  - EmergencyBraking/DesignTimeVerification/ (reach_exact.m, reach_approx.m, etc.)
  - Sherlock-Benchmarks/ (13 subdirectories, all use LayerS)
  - Buck_Hybrid/ (main_reach.m)
  - Tutorial/NNCS/AEBS/ (uses FFNN, LayerS)

#### 12. Add README Files to Examples
- **Priority**: MEDIUM
- **Effort**: 1-2 weeks
- **Issue**: Only 2 README files found in entire examples/ structure

**Directories needing READMEs**:
- NN/ACASXU/
- NN/MNIST/
- NN/NeuralODEs/
- NN/SemanticSegmentation/
- All NNCS subdirectories (except Tutorial)
- NN/ImperialCollege/
- NN/SHERLOCK/
- NN/cifar10/, medmnist/, xai/

---

## Medium Priority - Short Term (1-3 months)

### Missing Features - Modern Neural Network Support

#### 13. Add Attention Mechanism Layers
- **Priority**: HIGH (for modern architectures)
- **Effort**: 4-6 weeks
- **Layers Needed**:
  - MultiHeadAttentionLayer
  - SelfAttentionLayer
  - CrossAttentionLayer
  - PositionalEncodingLayer
- **Impact**: Critical for Transformer verification
- **Approach**: Start with over-approximate reachability

#### 14. Add Modern Normalization Layers
- **Priority**: MEDIUM-HIGH
- **Effort**: 2-3 weeks
- **Layers Needed**:
  - GroupNormalizationLayer (high priority - common in vision)
  - InstanceNormalizationLayer
  - RMSNormalizationLayer (used in modern transformers)

#### 15. Add Modern Activation Functions
- **Priority**: MEDIUM
- **Effort**: 2-3 weeks each
- **Functions Needed**:
  - GELU (critical for transformers)
  - Swish/SiLU
  - Mish
  - PReLU (parametric ReLU)

#### 16. Add Specialized Layers
- **Priority**: MEDIUM
- **Effort**: 1-2 weeks each
- **Layers**:
  - SeparableConv2DLayer (depthwise + pointwise)
  - DilatedConv2DLayer (atrous convolution)
  - SqueezeExcitationLayer (SE blocks)

### Incomplete Features

#### 17. Add Normalization Options to Input Layers
- **Priority**: MEDIUM
- **Effort**: 2-3 weeks
- **Files**:
  - `ImageInputLayer.m:83`
  - `Image3DInputLayer.m:75`
  - `FeatureInputLayer.m:93`
- **Issue**: "TODO: add normalization options"
- **Impact**: Limited preprocessing support

#### 18. Implement Feedback Support in NNCS
- **Priority**: MEDIUM
- **Effort**: 2-4 weeks
- **Files**:
  - `LinearNNCS.m:6`
  - `DLinearNNCS.m:6`
- **Issue**: "TODO: need support feedback y(k), y(k-d)"
- **Impact**: Limits control system verification capabilities

#### 19. Complete VolumeStar Helper Methods
- **Priority**: LOW-MEDIUM
- **Effort**: 1-2 weeks
- **File**: `VolumeStar.m`
- **Lines**: 888 ("check permuted array is correctly reshaped"), 935 ("begin here again")
- **Issue**: Indicates incomplete development

### Progress Tracking & Usability

#### 20. Add Progress Callbacks
- **Priority**: MEDIUM
- **Effort**: 1 week
- **Implementation**: Add callback support to NN.reach()

```matlab
% Add to NN.m
properties
    progressCallback = []
end

% In reach method:
if ~isempty(obj.progressCallback)
    obj.progressCallback(struct('layer', i, 'total', obj.numLayers, ...
                                 'name', obj.Layers{i-1}.Name));
end
```

#### 21. Add Verification Checkpointing
- **Priority**: MEDIUM
- **Effort**: 1 week
- **Purpose**: Save intermediate results to avoid losing progress on crashes
- **Implementation**: Add checkpointDir option to save layer-wise results

#### 22. Improve ONNX Import Logging
- **Priority**: MEDIUM
- **Effort**: 2-3 days
- **File**: `onnx2nnv.m`
- **Issue**: 7 nested try-catch blocks but no feedback on which method succeeded
- **Fix**: Add logging to show which import method worked and which ONNX ops were converted

---

## Long Term Improvements (3-12 months)

### Advanced Verification Features

#### 23. Implement AutoVerifier
- **Priority**: HIGH (strategic)
- **Effort**: 2-3 months
- **Purpose**: Automatically select best verification strategy based on network structure
- **Features**:
  - Analyze network complexity
  - Choose optimal reachability method
  - Timeout handling with fallback strategies
  - Adaptive exact/approximate switching

#### 24. Compositional Verification Framework
- **Priority**: HIGH (strategic)
- **Effort**: 3-6 months
- **Purpose**: Scale to very large networks via decomposition
- **Features**:
  - Network decomposition heuristics
  - Sub-network verification
  - Specification chaining
  - Result composition

#### 25. Complete GPU Acceleration
- **Priority**: MEDIUM-HIGH
- **Effort**: 2-3 months
- **Files**: Multiple "TODO: explore GPU" comments in:
  - `MaxPooling2DLayer.m:383`
  - `Conv2DLayer.m:916`
  - `Conv1DLayer.m:455`
  - `AveragePooling2DLayer.m:358`
- **Features**:
  - GPU-compatible LP solver
  - Batch LP solving on GPU
  - Layer reach methods with GPU paths

#### 26. Add Quantized Network Verification
- **Priority**: MEDIUM
- **Effort**: 3-4 months
- **Purpose**: Verify INT8/INT4 quantized networks
- **Impact**: Critical for edge deployment verification

#### 27. SMT Solver Integration
- **Priority**: MEDIUM
- **Effort**: 2-3 months
- **Purpose**: Complement LP-based methods with SMT solving
- **Impact**: Better handling of non-linear constraints

### Enhanced VNNLIB/ONNX Support

#### 28. Enhance VNNLIB Export
- **Priority**: MEDIUM
- **Effort**: 3-4 weeks
- **File**: `export2vnnlib.m`
- **Current Limitations**:
  - Only one input set allowed
  - Only one HalfSpace output property
  - Cannot export complex verification problems
- **Needed**: Support for disjunctions, multiple properties, non-linear constraints

#### 29. Add VNNLIB Roundtrip Testing
- **Priority**: LOW-MEDIUM
- **Effort**: 1 week
- **Purpose**: Ensure export/import consistency
- **Implementation**: Test suite that exports then re-imports properties

### Interactive Tools

#### 30. Create Interactive Visualization Tool
- **Priority**: LOW-MEDIUM
- **Effort**: 1-2 months
- **Features**:
  - 2D/3D reachable set visualization
  - Layer-by-layer animation
  - Interactive projections and slicing
  - Zoom and pan controls

#### 31. Add Verification Time Estimation
- **Priority**: MEDIUM
- **Effort**: 2-3 weeks
- **Purpose**: "Dry run" mode to estimate verification time
- **Implementation**: Quick conservative analysis without full verification

---

## Testing and Quality

### Write Missing Tests

#### 32. Complete todo_test.m Files (8 files)
- **Priority**: HIGH
- **Effort**: 2-3 weeks
- **Files**:
  - `tests/nn/layers/ODEblockLayer/todo_test.m`
  - `tests/nn/layers/TransposeConv1DLayer/todo_test.m`
  - `tests/nn/layers/Conv3DLayer/todo_test.m`
  - `tests/nn/layers/GlobalAveragePooling1DLayer/todo_test.m`
  - `tests/nn/layers/ElementwiseAffineLayer/todo_test.m`
  - `tests/nn/layers/Conv1DLayer/todo_test.m`
  - `tests/nn/layers/AveragePooling3DLayer/todo_test.m`
  - `tests/nn/layers/BatchNormalizationLayer/todo_test.m`

#### 33. Add Tests for Layers Without Test Files
- **Priority**: MEDIUM
- **Effort**: 1-2 weeks
- **Layers**:
  - ImageInputLayer
  - LayerS
  - MaxUnpooling2DLayer

#### 34. Add Integration Tests
- **Priority**: HIGH
- **Effort**: 2 weeks
- **Purpose**: End-to-end verification workflows
- **Tests Needed**:
  - MNIST verification workflow
  - VNNLIB roundtrip (export then import)
  - ONNX import → verify → export
  - NNCS safety verification

#### 35. Add Performance Regression Tests
- **Priority**: MEDIUM
- **Effort**: 1 week
- **Purpose**: Ensure optimizations don't regress
- **Implementation**: Benchmark suite with timing checks

#### 36. Add Fuzzing Tests
- **Priority**: LOW-MEDIUM
- **Effort**: 1-2 weeks
- **Purpose**: Random network/input testing for robustness
- **Implementation**: Generate random networks and verify no crashes

#### 37. Fix VolumeStar Test Performance Issues
- **Priority**: LOW
- **Effort**: 1-2 weeks
- **Files**: `test_VolumeStar.m`, `test_gpu_VolumeStar.m`
- **Issue**: Tests disabled due to:
  - Line 26: `sample(10)` takes too long
  - Line 41: `contains()` causes out of memory
- **Fix**: Optimize underlying methods

---

## Documentation

### Code Documentation

#### 38. Document Activation Function Methods
- **Priority**: LOW-MEDIUM
- **Effort**: 1-2 weeks
- **Files**: LeakyReLU.m, Sign.m, PosLin.m, etc.
- **Issue**: Most methods have minimal header comments
- **Needed**: Add examples, performance notes, algorithm references

#### 39. Add Inline Examples to Key Methods
- **Priority**: MEDIUM
- **Effort**: 1 week
- **Files**: NN.m (evaluate, reach, verify methods)
- **Needed**: Usage examples in method documentation

#### 40. Create API Reference Documentation
- **Priority**: MEDIUM
- **Effort**: 2-3 weeks
- **Implementation**: Auto-generate with MATLAB's doc system or Sphinx
- **Publish**: HTML documentation on GitHub Pages

#### 41. Create Architecture Guide
- **Priority**: LOW-MEDIUM
- **Effort**: 1 week
- **Content**: Expand CLAUDE.md into developer documentation
- **Include**: Class hierarchies, design patterns, algorithm explanations

#### 42. Add Troubleshooting Guide
- **Priority**: MEDIUM
- **Effort**: 1 week
- **Content**: Common errors and solutions
- **Format**: FAQ-style with search keywords

#### 43. Create Performance Tuning Guide
- **Priority**: MEDIUM
- **Effort**: 1 week
- **Content**: Guide for choosing methods, parameters, parallel options
- **Include**: Benchmarks, decision flowchart

#### 44. Add Contributor Guide
- **Priority**: LOW
- **Effort**: 3-4 days
- **Content**: Coding standards, PR process, testing requirements

---

## Performance Optimizations

### Critical Performance Issues (10-50x speedup potential)

#### 45. Implement Batch LP Solver Calls
- **Priority**: VERY HIGH
- **Effort**: 3-4 weeks
- **Impact**: Expected 10-50x speedup
- **Files**: Star.m, ImageStar.m, lpsolver.m
- **Issue**: Sequential LP solving in loops (getBox, getRanges, reach methods)
- **Implementation**: Combine multiple LP problems into single batch call

#### 46. Cache Bounds in Star and ImageStar
- **Priority**: VERY HIGH
- **Effort**: 1-2 weeks
- **Impact**: Expected 5-20x speedup for repeated queries
- **Files**: Star.m, ImageStar.m
- **Implementation**:
  - Check if state_lb/state_ub already computed before calling LP solver
  - Lazily compute and cache im_lb/im_ub in ImageStar
  - Invalidate cache when set is modified

#### 47. Fix Array Growing in Loops
- **Priority**: HIGH
- **Effort**: 1 week
- **Impact**: Expected 2-5x speedup
- **Files**:
  - `Star.m:296-305` (sample method)
  - `Star.m:1880-1883` (get_hypercube_hull)
  - `Zono.m:266-274` (generator sorting)
  - `PosLin.m:279-281, 285-287` (stepReachMultipleInputs)
- **Fix**: Pre-allocate arrays

### High Impact Optimizations (2-5x speedup)

#### 48. Vectorize Zonotope Operations
- **Priority**: HIGH
- **Effort**: 1 week
- **File**: `Zono.m:266-274`
- **Current**: Loop to compute generator norms
- **Fix**: Use `vecnorm(obj.V, 2, 1)`

#### 49. Minimize Type Conversions in lpsolver
- **Priority**: MEDIUM-HIGH
- **Effort**: 3-5 days
- **File**: `lpsolver.m:23-26, 38-41`
- **Issue**: Unnecessary GPU→CPU and single→double conversions
- **Fix**: Allow solvers to work with single precision; batch GPU transfers

#### 50. Default to Parallel for Large Problems
- **Priority**: MEDIUM
- **Effort**: 1 week
- **Files**: Star.m (getMins, getMaxs), layer reach methods
- **Implementation**: Automatically use parallel when dim > threshold (e.g., 100)

### Medium Impact Optimizations (1.5-3x speedup)

#### 51. Eliminate Unnecessary Set Conversions
- **Priority**: MEDIUM
- **Effort**: 2-3 weeks
- **Files**: ReluLayer.m and other activation layers
- **Issue**: Repeated Star↔ImageStar/VolumeStar conversions
- **Fix**: Implement reach methods directly on ImageStar/VolumeStar

#### 52. Implement GPU-Compatible LP Solver
- **Priority**: MEDIUM
- **Effort**: 1-2 months
- **Impact**: Expected 2-10x speedup for GPU workflows
- **Approach**: Integrate cuOpt or custom GPU LP implementation

#### 53. Optimize MaxPooling Exact Reachability
- **Priority**: MEDIUM
- **Effort**: 2-3 weeks
- **File**: `MaxPooling2DLayer.m:556-588`
- **Issue**: Split operations could be parallelized better
- **Fix**: Process independent splits in parallel

#### 54. Optimize Star.sample() Method
- **Priority**: LOW
- **Effort**: 1 week
- **File**: `Star.m:277`
- **Comment**: "TODO: optimize sampling"

#### 55. Minimize LP Calls in LeakyReLU
- **Priority**: LOW-MEDIUM
- **Effort**: 2-3 weeks
- **File**: `LeakyReLU.m:1746`
- **Comment**: "TODO: minimize the number of LP needs to be solved"

---

## Refactoring Opportunities

### High Impact Refactoring

#### 56. Extract Reach Options Parsing Helper
- **Priority**: HIGH
- **Effort**: 2-3 days
- **File**: NN.m
- **Lines**: 201-238 and 401-420 (duplicated code)
- **Fix**: Extract to `parseReachOptions(reachOptions)` helper method

#### 57. Create PoolingLayer Base Class
- **Priority**: HIGH
- **Effort**: 1 week
- **Files**: MaxPooling2DLayer.m, AveragePooling2DLayer.m
- **Issue**: ~150 lines of duplicated validation and getter/setter methods
- **Fix**: Extract shared logic to base class

#### 58. Refactor Star Constructor
- **Priority**: HIGH
- **Effort**: 2 weeks
- **File**: `Star.m:29-275` (246 lines)
- **Issue**: Massive switch statement, hard to maintain
- **Fix**: Extract to named constructors:
  - `Star.fromVCD(V, C, d, ...)`
  - `Star.fromBounds(lb, ub)`
  - `Star.fromPolyhedron(P)`

#### 59. Refactor ImageStar Constructor
- **Priority**: HIGH
- **Effort**: 2 weeks
- **File**: `ImageStar.m:99-365` (266 lines)
- **Issue**: Similar to Star constructor
- **Fix**: Extract to named constructors

#### 60. Simplify NN.verify_segmentation
- **Priority**: MEDIUM-HIGH
- **Effort**: 1 week
- **File**: `NN.m:835-991` (156 lines)
- **Issue**: Complex nested loops
- **Fix**: Extract helper methods for pixel classification and metrics

### Medium Impact Refactoring

#### 61. Create LP Solver Result Checking Utility
- **Priority**: MEDIUM
- **Effort**: 3-5 days
- **Files**: Star.m, ImageStar.m (100+ instances)
- **Issue**: Repeated exitflag checking pattern
- **Fix**: `checkLPResult(fval, exitflag, message)` helper

#### 62. Create processMultipleInputs Helper
- **Priority**: MEDIUM
- **Effort**: 1 week
- **Files**: MaxPooling2DLayer.m, AveragePooling2DLayer.m, etc.
- **Issue**: Identical parallel/single core branching logic
- **Fix**: Base class method `processMultipleInputs(inputs, method, option)`

#### 63. Standardize Method Naming
- **Priority**: MEDIUM (breaking change)
- **Effort**: 2-3 weeks
- **Issue**: Mixed snake_case and camelCase
- **Examples**:
  - `reach_star_single_input` → `reachStarSingleInput`
  - `get_localMax_index` → `getLocalMaxIndex`
  - `set_poolSize` → `setPoolSize`
- **Note**: Breaking change, requires deprecation warnings

#### 64. Simplify NN Constructor
- **Priority**: MEDIUM
- **Effort**: 1 week
- **File**: `NN.m:54-108`
- **Issue**: Deep nesting with 5 switch cases
- **Fix**: Use builder pattern or parameter struct

#### 65. Use MATLAB Arguments Block (Requires R2021a+)
- **Priority**: LOW-MEDIUM
- **Effort**: 2-3 weeks
- **Files**: Multiple reach methods with nargin switching
- **Fix**: Use modern arguments validation
- **Note**: Requires MATLAB version bump

```matlab
function IS = reach(obj, in_images, method, options)
    arguments
        obj
        in_images
        method = 'approx-star'
        options.parallel = []
        options.display = []
        options.lp_solver = 'linprog'
    end
end
```

### Code Cleanup

#### 66. Remove Commented-Out Code
- **Priority**: LOW
- **Effort**: 2-3 days
- **Files**:
  - `ImageStar.m:513-527` (entire recurrentMap method commented)
  - `Star.m:754-760` (progress tracking code)
  - Multiple "TODO: explore GPU" comments
- **Action**: Either implement or remove

#### 67. Define Magic Number Constants
- **Priority**: LOW-MEDIUM
- **Effort**: 1 week
- **Examples**:
  - Star.m line 366, 402: `0.0001` tolerance → `FLOAT_TOLERANCE` constant
  - LeakyReluLayer.m line 9: `gamma = 0.01` → named constant
  - NN.m defaults → `DEFAULT_REACH_METHOD`, `DEFAULT_RELAX_FACTOR`

#### 68. Remove LayerS and FFNNS Deprecation
- **Priority**: LOW (long-term)
- **Effort**: 1 week
- **Action**: After all examples updated, add deprecation warnings to old classes
- **Eventually**: Remove deprecated classes (major version bump)

---

## Installation & Build

#### 69. Make GLPK Patching More Robust
- **Priority**: MEDIUM
- **Effort**: 1 week
- **File**: `install.m:37-54`
- **Issue**: Fragile file editing breaks if GLPK version/path changes
- **Fix**: Provide pre-patched GLPK or use wrapper functions

#### 70. Add Version Checking to Installation
- **Priority**: MEDIUM
- **Effort**: 3-5 days
- **File**: install.m
- **Checks Needed**:
  - MATLAB version >= 2023a
  - Required toolboxes installed
  - Python version for probabilistic reachability
  - Submodule initialization status

#### 71. Add Python Environment Verification
- **Priority**: LOW-MEDIUM
- **Effort**: 2-3 days
- **File**: install_ubuntu.sh
- **Issue**: Creates .venv but no verification it works
- **Fix**: Test imports after environment creation

---

## Summary Statistics

### By Priority

| Priority | Count | Estimated Total Effort |
|----------|-------|------------------------|
| 🔴 CRITICAL | 1 | 5 minutes |
| 🟠 HIGH | 28 | 25-35 weeks |
| 🟡 MEDIUM | 33 | 30-50 weeks |
| 🟢 LOW | 9 | 10-15 weeks |
| **TOTAL** | **71** | **65-100 weeks** |

### By Category

| Category | Count |
|----------|-------|
| Bug Fixes | 3 |
| User Experience | 8 |
| Missing Features | 13 |
| Example Updates | 5 |
| Testing | 7 |
| Documentation | 7 |
| Performance | 11 |
| Refactoring | 13 |
| Installation | 4 |

### Near-Term Actionable Items (Next 4 weeks)

1. ✅ Fix GlobalAveragePooling1DLayer bug (5 min)
2. ✅ Add NN.parse() static method (1 day)
3. ✅ Create ReachOptions class (2-3 days)
4. ✅ Improve error messages (1 week)
5. ✅ Add install_check.m (2 days)
6. ✅ Update NNCS/ACC examples (2-3 days)
7. ✅ Update NN/MNIST examples (1 week)
8. ✅ Complete UpsampleLayer reach (3-5 days)

**Total Near-Term Effort**: ~3-4 weeks

---

## Notes

- **Breaking Changes**: Items marked as breaking changes (like method renaming) should be done together in a major version release
- **MATLAB Version**: Some improvements (like arguments block) require MATLAB R2021a+
- **Testing First**: For any refactoring, ensure comprehensive tests exist first
- **Documentation**: Update CLAUDE.md and docs when implementing features
- **Examples**: Keep tutorial examples in sync with API changes

---

**Contributing**: See CLAUDE.md for development conventions and git workflow.
**Issues**: Report bugs at https://github.com/verivital/nnv/issues
**Contact**: Diego Manzanas Lopez (diego.manzanas.lopez@vanderbilt.edu)
