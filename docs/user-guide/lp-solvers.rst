LP Solvers & Configuration
==========================

.. rst-class:: lead

   Star-based reachability relies on linear programming (LP) to compute bounds
   and check containment. NNV supports three LP solver backends.

----

Supported Solvers
-----------------

.. list-table::
   :header-rows: 1
   :widths: 20 25 25 30

   * - Solver
     - License
     - Speed
     - Notes
   * - ``linprog``
     - MATLAB built-in
     - Baseline
     - Always available; requires Optimization Toolbox
   * - ``glpk``
     - Open source (GPL)
     - Faster
     - Included with NNV; good default alternative
   * - ``gurobi``
     - Commercial (free academic)
     - Fastest
     - Requires separate installation and license

Selecting a Solver
------------------

Set the solver in the ``reachOptions`` structure:

.. code-block:: matlab

   reachOptions.lp_solver = 'linprog';  % MATLAB built-in
   reachOptions.lp_solver = 'glpk';     % GNU Linear Programming Kit
   reachOptions.lp_solver = 'gurobi';   % Gurobi (if installed)

The ``lpsolver()`` function provides a unified interface used internally
by all Star set operations. It handles:

- Equality and inequality constraints
- Solver-specific option formatting
- Automatic fallback between solvers
- GPU array support

Performance Considerations
--------------------------

- **Gurobi** is significantly faster for large-scale verification problems
  and is the solver of choice for VNN-COMP submissions
- **GLPK** offers a good balance of speed and availability, included with NNV
- **linprog** is the fallback and always available with the Optimization Toolbox
- For ``approx-zono`` and ``abs-dom`` methods, LP solving is not required, making
  them faster regardless of solver choice

Installing Gurobi
------------------

1. Download from `gurobi.com <https://www.gurobi.com/>`_ (free academic license)
2. Install and configure the MATLAB interface
3. Verify with:

.. code-block:: matlab

   gurobi_version = gurobi_mex('version');
