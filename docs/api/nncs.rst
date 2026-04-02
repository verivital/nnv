Control Systems
===============

.. rst-class:: lead

   Plant models and closed-loop neural network control system classes.

----

Plant Models
------------

LinearODE
^^^^^^^^^

Continuous-time linear system: x' = Ax + Bu, y = Cx + Du.

.. code-block:: matlab

   plant = LinearODE(A, B, C, D)
   plant = LinearODE(A, B, C, D, controlPeriod, numReachSteps)

**Key Properties:** ``A``, ``B``, ``C``, ``D``, ``dim``, ``nI``, ``nO``,
``controlPeriod``, ``numReachSteps``

DLinearODE
^^^^^^^^^^

Discrete-time linear system: x[k+1] = Ax[k] + Bu[k].

.. code-block:: matlab

   plant = DLinearODE(A, B, C, D)

NonLinearODE
^^^^^^^^^^^^

Continuous-time nonlinear system via CORA. x' = f(x, u, t), y = Cx.

.. code-block:: matlab

   plant = NonLinearODE(nVars, nInputs, @dynamics_func, reachStep, controlPeriod, outputMat)

**CORA Options:** ``taylorTerms``, ``zonotopeOrder``, ``maxError``, ``tensorOrder``,
``alg``, ``reductionTechnique``

DNonLinearODE
^^^^^^^^^^^^^

Discrete-time nonlinear system.

.. code-block:: matlab

   plant = DNonLinearODE(nVars, nInputs, @dynamics_func, outputMat)

HybridA
^^^^^^^

Hybrid automaton (switched systems) via CORA.

.. code-block:: matlab

   plant = HybridA(nModes, nVars, nInputs, dynamics, ...)

NNCS Composition
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 25 25 20

   * - Class
     - Plant Type
     - Time
     - Constructor
   * - ``LinearNNCS``
     - LinearODE
     - Continuous
     - ``LinearNNCS(controller, plant)``
   * - ``NonlinearNNCS``
     - NonLinearODE
     - Continuous
     - ``NonlinearNNCS(controller, plant)``
   * - ``DLinearNNCS``
     - DLinearODE
     - Discrete
     - ``DLinearNNCS(controller, plant)``
   * - ``DNonlinearNNCS``
     - DNonLinearODE
     - Discrete
     - ``DNonlinearNNCS(controller, plant)``
   * - ``HybridANNCS``
     - HybridA
     - Hybrid
     - ``HybridANNCS(controller, plant)``

**Common Properties:** ``feedbackMap``, ``ref_I`` (reference input set), ``init_set``

**Common Methods:**

.. code-block:: matlab

   [reachSets, simTraces, controlTraces] = nncs.reach(reachPRM)

where ``reachPRM`` contains: ``init_set``, ``ref_input``, ``numSteps``, ``reachMethod``, ``numCores``
