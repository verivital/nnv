# NNV Test Suite

**Full documentation:** [Testing & CI/CD](https://verivital.github.io/nnv/developer/testing.html)

## Quick Start

```matlab
cd("code/nnv");
install;
addpath(fullfile(pwd, 'tests', 'test_utils'));

% Run ALL tests with baseline comparison
[test_results, regression_results] = run_tests_with_regression('tests', 'compare', true, 'verbose', true);
```

This runs **833 tests** (full suite) or **470 tests** (quick mode).

See the [documentation](https://verivital.github.io/nnv/developer/testing.html) for test categories, baseline management, coverage metrics, CI/CD workflows, and troubleshooting.
