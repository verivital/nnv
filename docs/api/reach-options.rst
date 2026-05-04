reachOptions Reference
======================

.. rst-class:: lead

   Complete reference for the ``reachOptions`` struct passed to ``NN.reach()``,
   ``GNN.reach()``, and NNCS reachability methods.

----

All reachability options are passed as fields of a MATLAB struct:

.. code-block:: matlab

   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';
   reachOptions.numCores = 1;
   output_sets = net.reach(input_set, reachOptions);

Fields
------

.. list-table::
   :header-rows: 1
   :widths: 22 12 15 51

   * - Field
     - Type
     - Default
     - Description
   * - ``reachMethod``
     - string
     - ``'approx-star'``
     - Reachability method: ``'exact-star'``, ``'approx-star'``, ``'approx-zono'``, ``'abs-dom'``, ``'direct'``, ``'krylov'``, ``'ode45'``
   * - ``numCores``
     - int
     - ``1``
     - Number of parallel workers for computation
   * - ``relaxFactor``
     - double
     - ``0``
     - For approx-star: 0 = full LP (most precise), 1 = no LP (fastest)
   * - ``lp_solver``
     - string
     - ``'linprog'``
     - LP solver backend: ``'linprog'``, ``'glpk'``, ``'gurobi'``
   * - ``dis_opt``
     - string
     - ``[]``
     - Set to ``'display'`` for verbose debugging output
   * - ``device``
     - string
     - ``'cpu'``
     - Computation device: ``'cpu'`` or ``'gpu'``

Conformal Prediction Fields
----------------------------

Additional fields for probabilistic verification via ``verify_robustness_cp``:

.. list-table::
   :header-rows: 1
   :widths: 22 12 15 51

   * - Field
     - Type
     - Default
     - Description
   * - ``coverage``
     - double
     - ``0.99``
     - Coverage probability (1 - miscoverage level)
   * - ``confidence``
     - double
     - ``0.99``
     - Confidence level for the PAC guarantee
   * - ``train_mode``
     - string
     - ``'Linear'``
     - Surrogate model type: ``'Linear'`` or ``'ReLU'``
   * - ``train_device``
     - string
     - ``'cpu'``
     - Python training device: ``'cpu'`` or ``'gpu'``
   * - ``train_epochs``
     - int
     - ``100``
     - Number of surrogate training iterations
   * - ``train_lr``
     - double
     - ``0.001``
     - Surrogate learning rate

NNCS Reach Parameters
---------------------

For NNCS ``reach()`` methods, use a ``reachPRM`` struct:

.. list-table::
   :header-rows: 1
   :widths: 22 12 15 51

   * - Field
     - Type
     - Default
     - Description
   * - ``init_set``
     - Star
     - (required)
     - Initial state set
   * - ``ref_input``
     - Star
     - (required)
     - Reference input set (e.g., target setpoint)
   * - ``numSteps``
     - int
     - (required)
     - Number of control periods to simulate
   * - ``reachMethod``
     - string
     - ``'approx-star'``
     - Same options as NN reachability
   * - ``numCores``
     - int
     - ``1``
     - Parallel workers
