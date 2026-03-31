Neural Network Control Systems
===============================

.. rst-class:: lead

   NNV verifies closed-loop control systems where a neural network controller
   interacts with a physical plant model. This page covers all plant types,
   NNCS composition, and reachability workflows.

----

Overview
--------

A Neural Network Control System (NNCS) consists of:

1. **Plant**: The physical system dynamics (linear, nonlinear, or hybrid)
2. **Controller**: A neural network that computes control inputs from plant outputs
3. **Feedback loop**: The controller observes plant state and applies actions at
   discrete time steps

.. code-block:: text

   ┌──────────┐   u(t)   ┌───────┐   y(t)   ┌──────────────┐
   │    NN    │ ────────> │ Plant │ ────────> │   Feedback   │
   │ Controller│          │       │           │   Mapping    │
   └──────────┘          └───────┘           └──────┬───────┘
        ^                                           │
        │              observation                   │
        └───────────────────────────────────────────┘

Plant Models
------------

LinearODE
^^^^^^^^^

Continuous-time linear system :math:`\dot{x} = Ax + Bu`, :math:`y = Cx + Du`:

.. code-block:: matlab

   plant = LinearODE(A, B, C, D);
   plant.controlPeriod = 0.1;   % Control step size (seconds)
   plant.numReachSteps = 2;     % ODE steps per control period

DLinearODE
^^^^^^^^^^

Discrete-time linear system :math:`x_{k+1} = Ax_k + Bu_k`, :math:`y_k = Cx_k + Du_k`:

.. code-block:: matlab

   plant = DLinearODE(A, B, C, D);

NonLinearODE
^^^^^^^^^^^^

Continuous-time nonlinear system :math:`\dot{x} = f(x, u, t)`, :math:`y = Cx`.
Wraps CORA's ``nonlinearSys`` class for zonotope-based reachability:

.. code-block:: matlab

   plant = NonLinearODE(nVars, nInputs, @dynamics_func, reachStep, controlPeriod, outputMat);

   % Configure CORA reachability parameters
   plant.taylorTerms = 4;        % Taylor expansion order
   plant.zonotopeOrder = 20;     % Zonotope complexity
   plant.maxError = ones(nVars, 1) * 0.001;  % Maximum allowed error

**Dynamics function format:**

.. code-block:: matlab

   function dx = dynamics_func(x, u)
       % x: state vector, u: control input vector
       dx(1,1) = x(2);
       dx(2,1) = -9.81 * sin(x(1)) + u(1);
   end

DNonLinearODE
^^^^^^^^^^^^^

Discrete-time nonlinear system :math:`x_{k+1} = f(x_k, u_k)`:

.. code-block:: matlab

   plant = DNonLinearODE(nVars, nInputs, @dynamics_func, outputMat);

HybridA (Hybrid Automaton)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Switched systems with multiple modes, transitions via guards and jumps.
Wraps CORA's ``HybridAutomaton`` class:

.. code-block:: matlab

   plant = HybridA(nModes, nVars, nInputs, dynamics, ...);

NNCS Composition Classes
-------------------------

NNV provides dedicated classes for each plant-controller combination:

.. list-table::
   :header-rows: 1
   :widths: 30 25 25 20

   * - Class
     - Plant Type
     - Time
     - Dynamics
   * - ``LinearNNCS``
     - LinearODE
     - Continuous
     - Linear
   * - ``NonlinearNNCS``
     - NonLinearODE
     - Continuous
     - Nonlinear
   * - ``DLinearNNCS``
     - DLinearODE
     - Discrete
     - Linear
   * - ``DNonlinearNNCS``
     - DNonLinearODE
     - Discrete
     - Nonlinear
   * - ``HybridANNCS``
     - HybridA
     - Continuous
     - Hybrid/Switched

Verification Workflow
---------------------

**Step 1: Define the neural network controller**

.. code-block:: matlab

   layers = {FullyConnectedLayer(W1, b1), ReLULayer(), FullyConnectedLayer(W2, b2)};
   controller = NN(layers);

**Step 2: Define the plant model**

.. code-block:: matlab

   % Example: 6-state ACC system (lead car + ego car)
   A = [...];  % 6x6 state matrix
   B = [...];  % 6x1 input matrix
   C = [...];  % 3x6 output matrix (relative distance, relative velocity, ego velocity)
   D = [...];  % 3x1 feedthrough
   plant = LinearODE(A, B, C, D);
   plant.controlPeriod = 0.1;
   plant.numReachSteps = 2;

**Step 3: Compose the NNCS**

.. code-block:: matlab

   nncs = LinearNNCS(controller, plant);
   nncs.feedbackMap = [0];  % Output feedback mapping

**Step 4: Set initial conditions and run reachability**

.. code-block:: matlab

   % Initial state region
   init_set = Star(lb_state, ub_state);

   % Reference input (e.g., desired following distance)
   ref_input = Star(lb_ref, ub_ref);

   % Reachability parameters
   reachPRM.init_set = init_set;
   reachPRM.ref_input = ref_input;
   reachPRM.numSteps = 50;           % Number of control steps
   reachPRM.reachMethod = 'approx-star';
   reachPRM.numCores = 1;

   [reachSets, simTraces, controlTraces] = nncs.reach(reachPRM);

**Step 5: Check safety**

.. code-block:: matlab

   % Safety specification: safe distance > 10 + 1.4 * ego_velocity
   % (Define as HalfSpace constraints on the output)
   safe = true;
   for t = 1:length(reachSets)
       for s = 1:length(reachSets{t})
           if ~verify_specification(reachSets{t}(s), safety_spec)
               safe = false;
               break;
           end
       end
   end

Key Properties
--------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Property
     - Description
   * - ``feedbackMap``
     - Maps plant outputs to controller inputs (index array)
   * - ``ref_I``
     - Reference input set (e.g., desired setpoint)
   * - ``init_set``
     - Initial state region for reachability
   * - ``controlPeriod``
     - Time between control actions (seconds)
   * - ``reachStep``
     - ODE integration step within each control period
   * - ``reachSetTree``
     - Hierarchical structure storing temporal reachable sets

CORA Integration
----------------

NNV's nonlinear plant models wrap the **CORA** toolbox
(COntinuous Reachability Analyzer) for computing reachable sets of ODEs:

- **Taylor expansion** methods compute over-approximations of nonlinear dynamics
- **Zonotope propagation** tracks the reachable set through time
- Configuration options: ``taylorTerms``, ``zonotopeOrder``, ``maxError``,
  ``tensorOrder``, ``reductionTechnique``

CORA is included as a submodule and initialized automatically when NNV is installed.

Example Applications
--------------------

- :doc:`/examples/nncs` -- ACC, AEBS, Inverted Pendulum, DC-DC Buck converter
- `Tutorial/NNCS/ACC/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NNCS/ACC>`_ -- Adaptive Cruise Control
- `Tutorial/NNCS/AEBS/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NNCS/AEBS>`_ -- Automated Emergency Braking
