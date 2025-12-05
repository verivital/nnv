# NNV Test Coverage Report

**Generated**: December 1, 2025
**Mode**: Quick (soundness + regression tests)
**NNV Version**: 3.0.0

---

## Summary

| Metric | Value |
|--------|-------|
| **Tests Run** | 470 |
| **Tests Passed** | 470 |
| **Tests Failed** | 0 |
| **Pass Rate** | 100% |
| **Source Files Tracked** | 148 |
| **Files Executed** | 72 |
| **Code Coverage** | 48.6% |
| **Runtime** | 10.2 minutes |

---

## Coverage by Directory

| Directory | Total Files | Executed | Coverage |
|-----------|-------------|----------|----------|
| `engine/nn/` | 65 | 56 | **86.2%** |
| `engine/set/` | 9 | 8 | **88.9%** |
| `engine/utils/` | 26 | 6 | 23.1% |
| `engine/nncs/` | 11 | 2 | 18.2% |
| `engine/hyst/` | 36 | 0 | 0.0% |
| **TOTAL** | **148** | **72** | **48.6%** |

### Coverage Analysis

**Well-Covered (>80%)**:
- Neural network layers (`nn/`) - Core verification functionality
- Set representations (`set/`) - Star, Box, Zono, ImageStar

**Moderate Coverage (20-50%)**:
- Utility functions (`utils/`) - Helper functions, some specialized

**Low Coverage (<20%)**:
- NNCS control systems - Specialized hybrid system verification
- Hyst/SpaceEx - External tool integration (0% - requires Java)

---

## Test Categories

### Soundness Tests (Layer Verification)
All layers verified for mathematical soundness:

| Layer Type | Tests | Status |
|------------|-------|--------|
| AdditionLayer | 11 | PASS |
| RMSNormLayer | 12 | PASS |
| SiLULayer | 12 | PASS |
| SwiGLULayer | 10 | PASS |
| ReluLayer | 6 | PASS |
| SigmoidLayer | 2 | PASS |
| TanhLayer | 2 | PASS |
| LeakyReluLayer | 5 | PASS |
| Conv2DLayer | 4 | PASS |
| Conv3DLayer | 4 | PASS |
| MaxPooling2DLayer | 5 | PASS |
| AveragePooling2DLayer | 6 | PASS |
| FullyConnectedLayer | 4 | PASS |
| FlattenLayer | 4 | PASS |
| BatchNormalizationLayer | 4 | PASS |
| ... and 40+ more layers | | PASS |

### Regression Tests
| Test Suite | Tests | Status |
|------------|-------|--------|
| Layer outputs | 10 | PASS |
| Reachability methods | 9 | PASS |
| ACAS-Xu VNN-LIB | 11 | PASS |
| MNIST FC | 8 | PASS |
| Neural ODE | 8 | PASS |
| NN basic | 9 | PASS |
| NNCS basic | 8 | PASS |
| RNN sequence | 8 | PASS |
| Segmentation | 8 | PASS |
| Verify robustness | 8 | PASS |
| Verify safety | 8 | PASS |
| Falsification | 8 | PASS |

### Set Operation Tests
| Set Type | Tests | Status |
|----------|-------|--------|
| Star | Multiple | PASS |
| Box | 7 | PASS |
| Zono | 6 | PASS |
| ImageStar | Multiple | PASS |
| ImageZono | 6 | PASS |
| VolumeStar | 5 | PASS |
| HalfSpace | 6 | PASS |

---

## Uncovered Files

The following files were not executed during testing (first 20):

1. `engine/adjust_glpk.m` - Install-time utility
2. `engine/hyst/src/matlab/SpaceExToStateflow.m` - Hybrid system tool
3. `engine/hyst/src/matlab/createConfigFromSpaceEx.m`
4. `engine/hyst/src/matlab/nonsemanticTranslation.m`
5. `engine/hyst/src/matlab/semanticTranslation.m`
6. `engine/hyst/src/matlab/simulationLoop.m`
7. `engine/hyst/src/matlab/translateAutomaton.m`
8. ... and 69 more (mostly Hyst/SpaceEx integration)

**Note**: Most uncovered files are in the `hyst/` directory which requires Java runtime and is used for hybrid automata translation from SpaceEx format.

---

## Solver Availability

| Solver | Available | Status |
|--------|-----------|--------|
| linprog (Optimization Toolbox) | Yes | OK |
| GLPK | Yes | OK |
| Gurobi | Yes | OK |

---

## How to Run Coverage Analysis

```matlab
% Quick mode (soundness + regression only)
results = track_coverage('quick');

% Full mode (all tests)
results = track_coverage('full');

% Generate HTML report
results = track_coverage('quick', 'report', true);

% View summary
fprintf('%s\n', results.summary);
```

---

## Recommendations

1. **Increase NNCS coverage**: Add more control system verification tests
2. **Add Hyst tests**: Create tests that don't require Java for basic functionality
3. **Utils coverage**: Add unit tests for utility functions
4. **Tutorial validation**: Ensure all tutorials run without errors

---

## Version History

| Date | Coverage | Tests | Notes |
|------|----------|-------|-------|
| 2025-12-01 | 48.6% | 470 | Initial NNV 3.0 coverage report |
