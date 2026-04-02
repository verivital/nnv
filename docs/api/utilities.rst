Utilities
=========

.. rst-class:: lead

   Network loading, format conversion, verification helpers, and diagnostic
   functions.

----

Network Loading & Conversion
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Function
     - Description
   * - ``net = matlab2nnv(Mnetwork)``
     - Convert MATLAB SeriesNetwork, DAGNetwork, dlnetwork, or LayerGraph to NNV NN object
   * - ``net = onnx2nnv(onnxFile)``
     - Import ONNX network directly to NNV format
   * - ``net = load_NN_from_mat(matFile)``
     - Load network from .mat file

Specification & Verification
-----------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Description
   * - ``[lb, ub, prop] = load_vnnlib(file)``
     - Parse VNNLIB file. Returns input bounds and output HalfSpace properties.
   * - ``result = verify_specification(reachSet, property)``
     - Check if reachable set satisfies property. Returns: 0=property failed (unsafe),
       1=property satisfied (safe), 2=unknown.
   * - ``export2vnnlib(lb, ub, prop, file)``
     - Export specification to VNNLIB format.

LP Solving
----------

.. code-block:: matlab

   [fval, exitflag] = lpsolver(f, A, b, Aeq, beq, lb, ub, solver, opts)

Unified LP interface supporting ``'linprog'``, ``'glpk'``, and ``'gurobi'``.
Handles automatic fallback between solvers.

Diagnostics & Info
------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Function
     - Description
   * - ``check_nnv_setup()``
     - Diagnostic tool: checks MATLAB version, toolboxes, submodules, Python env
   * - ``v = NNVVERSION()``
     - Returns version string (e.g., ``'NNV v3.0.0'``)
   * - ``path = nnvroot()``
     - Returns NNV installation root path
   * - ``python_path = cp_env()``
     - Returns path to Python executable for conformal prediction

Weight Perturbation (WPutils)
-----------------------------

``WPutils`` is a static class. All methods are called as ``WPutils.method_name(...)``.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Method
     - Description
   * - ``WPutils.get_weights_range(net, layer_no)``
     - Get the weight range (max - min) for a layer
   * - ``WPutils.print_layers_info(net)``
     - Display layer information for all layers
   * - ``WPutils.sample_weight_perturbed_nns(net, args)``
     - Sample networks with randomly perturbed weights

Weight perturbations are specified per-layer via the ``weightPerturb`` property
on ``FullyConnectedLayer`` and ``Conv2DLayer``. Each row of ``weightPerturb`` is
``[linear_index, lower_bound, upper_bound]``.

Probabilistic Verification
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Function
     - Description
   * - ``result = verify_robustness_cp(net, IS, reachOptions, target, nClasses)``
     - Conformal prediction robustness verification
   * - ``spec = CP_specification(coverage, confidence)``
     - Compute required sample sizes for CP verification
   * - ``reachSet = Prob_reach(net, IS, reachOptions)``
     - Main probabilistic reachability driver
