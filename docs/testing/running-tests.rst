Running Tests
=============

.. rst-class:: lead

   How to run NNV's test suite, manage baselines, and troubleshoot
   common issues.

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
