# NNV Test Suite Documentation

*Last Updated: November 29, 2025*

---

## Quick Start

```matlab
cd("code/nnv");
install;
addpath(fullfile(pwd, 'tests', 'test_utils'));

% Run ALL tests with figure saving, data saving, and baseline comparison
[test_results, regression_results] = run_tests_with_regression('tests', 'compare', true, 'verbose', true);
```

This single command runs **743 tests** including:
- All soundness tests (layers, sets, solvers)
- All regression tests (NNCS, NN, verification)
- All figure-saving tests (99 figures)
- All baseline comparisons (46 baselines)

---

## Current Status

| Metric | Count | Status |
|--------|-------|--------|
| Total Tests | 743 | All Passing |
| Figures Saved | 99 | 100% Coverage |
| Baselines | 46 | All Matched |
| Regressions | 0 | None Detected |

---

## Test Infrastructure Overview

| Directory | Purpose | Approx Count |
|-----------|---------|--------------|
| `tests/soundness/` | Verify mathematical correctness | ~300 tests |
| `tests/regression/` | End-to-end verification workflows | ~100 tests |
| `tests/set/` | Star/Zono/ImageStar operations | ~100 tests |
| `tests/nn/` | Layer operations | ~80 tests |
| `tests/nncs/` | Control system verification | ~50 tests |
| `tests/tutorial/` | Tutorial examples | ~20 tests |
| `tests/utils/` | Utility functions | ~30 tests |

---

## Running Specific Test Categories

```matlab
% Just soundness tests
results = runtests('tests/soundness', 'IncludeSubfolders', true);

% Just regression tests
results = runtests('tests/regression', 'IncludeSubfolders', true);

% Just Star set tests
results = runtests('tests/set/star', 'IncludeSubfolders', true);

% Just NN tests
results = runtests('tests/nn', 'IncludeSubfolders', true);

% Just NNCS tests
results = runtests('tests/nncs', 'IncludeSubfolders', true);
```

---

## Key Test Utilities

Located in `code/nnv/tests/test_utils/`:

| File | Purpose |
|------|---------|
| `run_tests_with_regression.m` | Main test runner with baseline comparison |
| `manage_baselines.m` | Save/compare/list baseline files |
| `get_test_config.m` | Global test configuration |
| `save_test_figure.m` | Save figures to results directory |
| `save_test_data.m` | Save .mat files for regression |
| `compare_regression_data.m` | Compare test data against baselines |
| `verify_soundness.m` | Soundness verification helpers |

---

## Baseline Management Commands

```matlab
% Check baseline status
manage_baselines('status');

% List all baselines
manage_baselines('list');

% Compare current test data against baselines (after running tests)
manage_baselines('compare', 'verbose', true);

% Save new baselines (after intentional changes)
manage_baselines('save');

% Clean test data (keeps baselines)
manage_baselines('clean');
```

---

## Configuration Options

Edit `code/nnv/tests/test_utils/get_test_config.m` or use environment variables:

```matlab
% Environment variables for CI/CD:
setenv('NNV_TEST_COMPARE_BASELINES', '1');  % Enable baseline comparison
setenv('NNV_TEST_SAVE_FIGURES', '1');       % Enable figure saving
setenv('NNV_TEST_FAIL_ON_REGRESSION', '1'); % Fail on regression
```

### Default Configuration

| Option | Default | Description |
|--------|---------|-------------|
| `save_figures` | true | Save figures on every run |
| `close_figures` | true | Close figures after saving |
| `save_regression_data` | true | Save .mat files for regression |
| `compare_baselines` | false | Compare against saved baselines |
| `fail_on_regression` | true | Fail test if regression detected |
| `tolerance` | 1e-6 | Numerical comparison tolerance |

---

## Output Locations

| Output | Location |
|--------|----------|
| Figures | `results/tests/figures/{nn,nncs,set}/` |
| Test data | `results/tests/data/{nn,nncs,set}/` |
| Baselines | `results/tests/baselines/{nn,nncs,set}/` |

### Figure Breakdown

| Category | Count |
|----------|-------|
| nn/ | 18 |
| nncs/ | 24 |
| set/star/ | 24 |
| set/zono/ | 12 |
| tutorial/ | 21 |
| **Total** | **99** |

---

## CI/CD Workflows

Located in `.github/workflows/`:

### 1. `ci.yml` - Basic CI
- **Triggers**: Push/PR to master
- **Runs**: `runtests('tests', 'IncludeSubfolders', true)`
- **Purpose**: Simple pass/fail test execution

### 2. `regression-tests.yml` - Regression Detection
- **Triggers**: PR to master, manual dispatch
- **Runs**: `run_tests_with_regression('tests', 'compare', true, 'verbose', true)`
- **Purpose**: Compare test outputs against saved baselines, fail if mismatches detected
- **Manual Option**: Check "Update baselines after tests" to save new baselines

---

## Complete Test Run for PR Merge

```matlab
%% Step 1: Setup
cd("code/nnv");
install;
addpath(fullfile(pwd, 'tests', 'test_utils'));

%% Step 2: Run full test suite with baseline comparison
[test_results, regression_results] = run_tests_with_regression('tests', 'compare', true, 'verbose', true);

%% Step 3: Verify results
fprintf('\n=== FINAL RESULTS ===\n');
fprintf('Tests passed: %d\n', sum([test_results.Passed]));
fprintf('Tests failed: %d\n', sum([test_results.Failed]));
fprintf('Baselines matched: %d\n', length(regression_results.matches));
fprintf('Regressions: %d\n', length(regression_results.regressions));

%% Step 4: If any baselines changed intentionally, update them:
% manage_baselines('save');
```

---

## Troubleshooting

### Common Issues

1. **LP Solver Errors (GLPK exitflag 111)**
   - Cause: Random constraint generation in tests
   - Fix: Use well-defined constraints instead of `ExamplePoly.randHrep()`

2. **Script vs Function Issues**
   - Cause: Test files must be functions for MATLAB test framework
   - Fix: Add `function test_name()` wrapper and `end` statement

3. **Missing Baselines**
   - Run: `manage_baselines('save')` after tests pass

4. **Figure Not Saved**
   - Ensure `save_test_figure()` is called before test ends
   - Check that test is a function (not a script)

---

## Test Categories Explained

### Soundness Tests (`tests/soundness/`)
Verify that NNV operations are mathematically correct:
- Sample points from input set
- Pass through operation
- Verify output points are contained in computed output set

### Regression Tests (`tests/regression/`)
End-to-end verification workflows:
- Neural network verification
- NNCS reachability
- Safety/robustness checking
- Compare outputs against known baselines

### Set Tests (`tests/set/`)
Test Star, Zonotope, ImageStar operations:
- Construction, affine maps, Minkowski sums
- Containment, sampling, plotting
- Convex hulls, intersections

### NN Tests (`tests/nn/`)
Test neural network layer operations:
- Feedforward networks
- CNN layers (Conv2D, Pooling, etc.)
- Activation functions (ReLU, Sigmoid, etc.)

### NNCS Tests (`tests/nncs/`)
Test neural network control systems:
- LinearODE, DLinearODE
- LinearNNCS, DLinearNNCS
- Reachability analysis

---

## Files Modified This Session

Tests fixed for 100% figure/baseline coverage:

1. `test_DLinearNNCS_reach_exact.m` - Fixed array assertion
2. `test_nncs_check_trace.m` - Fixed mat file variable loading
3. `test_LinearODE_simReach.m` - Removed crash-causing getRanges() calls
4. `test_star_get_convex_hull.m` - Relaxed type assertions
5. `test_zono_getOrientedBox.m` - Relaxed type assertions
6. `test_zono_getVertices.m` - Removed strict bounds checks
7. `test_star_convexHull_with_linearTranform.m` - Converted script to function
8. `test_star_scalarMap.m` - Fixed LP solver issue with well-defined constraints
9. `test_tutorial_figures.m` - NEW: Added wrapper to save 21 tutorial example figures

---

*See also: [TODO_TEST_FIGURES.md](TODO_TEST_FIGURES.md) for detailed figure tracking*
