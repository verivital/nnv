CI/CD Pipeline
==============

.. rst-class:: lead

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
