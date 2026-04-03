Testing & CI/CD
================

.. rst-class:: lead

   How to run NNV's test suite, manage baselines, check coverage,
   and work with the CI/CD pipeline.

----

Quick Start
-----------

.. code-block:: matlab

   cd("code/nnv");
   install;
   addpath(fullfile(pwd, 'tests', 'test_utils'));

   % Run ALL tests with regression comparison
   [test_results, regression_results] = run_tests_with_regression( ...
       'tests', 'compare', true, 'verbose', true);

This runs **833 tests** (full suite) or **470 tests** (quick mode) covering:

- Soundness tests (layers, sets, solvers)
- Regression tests (end-to-end verification workflows)
- Figure-saving tests (99 figures)
- Baseline comparisons (46 baselines)

Test Categories
---------------

.. list-table::
   :header-rows: 1
   :widths: 25 40 15

   * - Directory
     - Purpose
     - Count
   * - ``tests/soundness/``
     - Mathematical correctness verification
     - ~300
   * - ``tests/regression/``
     - End-to-end verification workflows
     - ~100
   * - ``tests/set/``
     - Star/Zono/ImageStar operations
     - ~100
   * - ``tests/nn/``
     - Layer operations
     - ~80
   * - ``tests/nncs/``
     - Control system verification
     - ~50
   * - ``tests/tutorial/``
     - Tutorial examples
     - ~20
   * - ``tests/utils/``
     - Utility functions
     - ~30

Running Specific Categories
---------------------------

.. code-block:: matlab

   % Just soundness tests
   results = runtests('tests/soundness', 'IncludeSubfolders', true);

   % Just Star set tests
   results = runtests('tests/set/star', 'IncludeSubfolders', true);

   % Just NN layer tests
   results = runtests('tests/nn', 'IncludeSubfolders', true);

   % Just NNCS tests
   results = runtests('tests/nncs', 'IncludeSubfolders', true);

Baseline Management
-------------------

.. code-block:: matlab

   % Check baseline status
   manage_baselines('status');

   % List all baselines
   manage_baselines('list');

   % Compare current test data against baselines
   manage_baselines('compare', 'verbose', true);

   % Save new baselines (after intentional changes)
   manage_baselines('save');

   % Clean test data (keeps baselines)
   manage_baselines('clean');

Configuration
-------------

Edit ``code/nnv/tests/test_utils/get_test_config.m`` or use environment variables:

.. code-block:: matlab

   setenv('NNV_TEST_COMPARE_BASELINES', '1');
   setenv('NNV_TEST_SAVE_FIGURES', '1');
   setenv('NNV_TEST_FAIL_ON_REGRESSION', '1');

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Default
     - Description
   * - ``save_figures``
     - true
     - Save figures on every run
   * - ``close_figures``
     - true
     - Close figures after saving
   * - ``save_regression_data``
     - true
     - Save .mat files for regression comparison
   * - ``compare_baselines``
     - false
     - Compare against saved baselines
   * - ``fail_on_regression``
     - true
     - Fail test if regression detected
   * - ``tolerance``
     - 1e-6
     - Numerical comparison tolerance

Output Locations
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Output
     - Location
   * - Figures
     - ``results/tests/figures/{nn,nncs,set}/``
   * - Test data
     - ``results/tests/data/{nn,nncs,set}/``
   * - Baselines
     - ``results/tests/baselines/{nn,nncs,set}/``

Figure Breakdown
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 20

   * - Category
     - Count
   * - nn/
     - 18
   * - nncs/
     - 24
   * - set/star/
     - 24
   * - set/zono/
     - 12
   * - tutorial/
     - 21
   * - **Total**
     - **99**

Test Categories Explained
--------------------------

**Soundness tests** (``tests/soundness/``):
Verify that NNV operations are mathematically correct -- sample points from
input set, pass through operation, verify output points are contained in
computed output set.

**Regression tests** (``tests/regression/``):
End-to-end verification workflows -- neural network verification, NNCS
reachability, safety/robustness checking. Compare outputs against known baselines.

**Set tests** (``tests/set/``):
Test Star, Zonotope, ImageStar operations -- construction, affine maps,
Minkowski sums, containment, sampling, plotting, convex hulls, intersections.

**NN tests** (``tests/nn/``):
Test neural network layer operations -- feedforward, CNN layers (Conv2D, Pooling),
activation functions (ReLU, Sigmoid).

**NNCS tests** (``tests/nncs/``):
Test neural network control systems -- LinearODE, DLinearODE, LinearNNCS,
DLinearNNCS, reachability analysis.

Troubleshooting
---------------

**LP Solver Errors (GLPK exitflag 111):**
Usually caused by random constraint generation. Use well-defined constraints
instead of ``ExamplePoly.randHrep()``.

**Script vs Function Issues:**
Test files must be MATLAB functions (not scripts) for the test framework.
Add ``function test_name()`` wrapper and ``end`` statement.

**Missing Baselines:**
Run ``manage_baselines('save')`` after a clean test pass.

**Figure Not Saved:**
Ensure ``save_test_figure()`` is called before the test ends.
Check that the test is a function (not a script).

Environment Variables
---------------------

.. code-block:: matlab

   setenv('NNV_TEST_COMPARE_BASELINES', '1');   % Enable baseline comparison
   setenv('NNV_TEST_SAVE_FIGURES', '1');         % Enable figure saving
   setenv('NNV_TEST_FAIL_ON_REGRESSION', '1');   % Fail on regression

Code Coverage
-------------

NNV 3.0 tracks code coverage via the ``track_coverage.m`` utility.

----

Current Coverage
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 20

   * - Directory
     - Coverage
   * - **Overall**
     - **48.6%**
   * - engine/nn/
     - 86.2%
   * - engine/set/
     - 88.9%
   * - engine/utils/
     - 23.1%
   * - engine/nncs/
     - 18.2%
   * - engine/hyst/
     - 0.0%

Running Coverage Analysis
-------------------------

.. code-block:: matlab

   addpath(fullfile(pwd, 'tests', 'test_utils'));
   coverage = track_coverage('tests', 'quick_mode', true);

The highest coverage is in the core verification components (``nn/`` and ``set/``),
which are the most critical for correctness.

See `TEST_COVERAGE.md <https://github.com/verivital/nnv/blob/master/TEST_COVERAGE.md>`_
for the detailed per-file coverage report.


   NNV uses GitHub Actions for automated testing on every push and
   pull request.

----

Workflows
---------

Located in ``.github/workflows/``:

ci.yml -- Basic CI
^^^^^^^^^^^^^^^^^^

- **Triggers**: Push/PR to master
- **Runs**: ``runtests('tests', 'IncludeSubfolders', true)``
- **Purpose**: Simple pass/fail test execution

regression-tests.yml -- Regression Detection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **Triggers**: PR to master, manual dispatch
- **Runs**: ``run_tests_with_regression('tests', 'compare', true, 'verbose', true)``
- **Purpose**: Compare test outputs against saved baselines, fail on mismatches
- **Manual option**: Check "Update baselines after tests" to save new baselines

Complete PR Merge Checklist
---------------------------

.. code-block:: matlab

   %% Step 1: Setup
   cd("code/nnv");
   install;
   addpath(fullfile(pwd, 'tests', 'test_utils'));

   %% Step 2: Run full test suite with baseline comparison
   [test_results, regression_results] = run_tests_with_regression( ...
       'tests', 'compare', true, 'verbose', true);

   %% Step 3: Verify results
   fprintf('\n=== FINAL RESULTS ===\n');
   fprintf('Tests passed: %d\n', sum([test_results.Passed]));
   fprintf('Tests failed: %d\n', sum([test_results.Failed]));
   fprintf('Baselines matched: %d\n', length(regression_results.matches));
   fprintf('Regressions: %d\n', length(regression_results.regressions));

   %% Step 4: If baselines changed intentionally:
   % manage_baselines('save');
