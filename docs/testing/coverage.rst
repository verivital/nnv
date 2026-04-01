Test Coverage
=============

.. rst-class:: lead

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
