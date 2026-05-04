Neural ODEs
===========

.. rst-class:: lead

   Verify continuous-time learned dynamics using Neural ODE layers.
   This tutorial covers defining ODE dynamics, configuring CORA
   integration parameters, and computing reachable state trajectories.

----

What You Will Learn
-------------------

- How Neural ODEs represent learned dynamics as continuous-time differential equations
- How to define a dynamics function and wrap it in an ``ODEblockLayer``
- How to configure CORA integration parameters (Taylor terms, zonotope order)
- How to compute and visualize reachable state trajectories over time
- When to use ``direct``, ``ode45``, or ``krylov`` reachability methods

Background
----------

A Neural ODE replaces discrete network layers with a continuous-time ODE:
instead of ``y = f(x)`` at each layer, the network solves ``dx/dt = f(x, t)``
over a time interval. This is useful for modeling physical dynamics,
continuous control policies, and time-series prediction.

NNV verifies Neural ODEs by computing the set of all possible state
trajectories starting from an uncertain initial condition, using the
CORA toolbox for ODE reachability.

Fixed-Point Attractor (FPA) Example
-------------------------------------

This example verifies a Continuous-Time Recurrent Neural Network that
converges to a fixed point attractor.

Prerequisites
^^^^^^^^^^^^^

- ``CTRNN_FPA.m`` -- dynamics function (included in the example folder)
- CORA toolbox (included as NNV submodule, initialized by ``install.m``)

Step 1: Define the Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dynamics function defines the ODE right-hand side. It must accept a
state vector and return the derivative:

.. code-block:: matlab

   function dx = CTRNN_FPA(x, u)
       % Continuous-Time Recurrent Neural Network dynamics
       % x: 5-dimensional state vector
       % u: input (unused in this example)
       % dx: state derivative
       %
       % Implements: dx/dt = -x + W*sigma(x) + b
       % where sigma is a sigmoid activation and W, b are learned parameters

       % ... (network-specific computation)
       dx = zeros(5, 1);
       dx(1) = ...;  % State derivatives
       % ...
   end

Step 2: Create the Neural ODE Layer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Define the nonlinear ODE system
   nVars = 5;            % State dimension
   nInputs = 1;          % Input dimension (can be 0 for autonomous systems)
   reachStep = 0.01;     % ODE integration step (seconds)
   controlPeriod = 10;   % Total simulation time (seconds)
   outputMat = eye(5);   % Output = full state (identity matrix)

   % Create NonLinearODE (wraps CORA's nonlinearSys)
   model = NonLinearODE(nVars, nInputs, @CTRNN_FPA, ...
       reachStep, controlPeriod, outputMat);

   % Configure CORA reachability parameters
   model.options.timeStep = 0.05;         % Verification time step (can be coarser than reachStep)
   model.options.taylorTerms = 4;         % Taylor expansion order (higher = more precise)
   model.options.zonotopeOrder = 20;      % Max zonotope generators (higher = more precise)
   model.options.alg = 'lin';             % Reachability algorithm
   model.options.tensorOrder = 2;         % Tensor approximation order

   % Wrap in an ODEblockLayer for use in an NN
   odelayer = ODEblockLayer(model);

.. admonition:: CORA Parameters Explained

   - **taylorTerms**: Order of the Taylor expansion used to approximate the
     nonlinear dynamics. Higher values reduce approximation error but
     increase computation. Typical: 4-6.
   - **zonotopeOrder**: Maximum number of generators in the zonotope
     representation. Controls the trade-off between precision and memory.
     Typical: 20-50.
   - **timeStep**: The verification step size. Can be larger than the ODE
     integration step for efficiency.
   - **alg**: Reachability algorithm. ``'lin'`` uses linearization-based
     approximation; other options available in CORA.

Step 3: Build the Neural ODE Network
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Create NNV network with ODE layer
   neuralode = NN({odelayer});

Step 4: Define Initial Set and Compute Reachability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Initial state with small uncertainty
   x0 = [0.21535; -0.58587; 0.8; 0.52323; 0.5];
   delta = 0.01;   % ±0.01 uncertainty in each state

   lb = x0 - delta;
   ub = x0 + delta;
   init_set = Star(lb, ub);

   % Compute reachable set over the full time horizon
   tic;
   R = neuralode.reach(init_set);
   fprintf('Neural ODE reachability completed in %.1f seconds\n', toc);

   % R is a cell array: R{t} contains the Star set(s) at time step t

Step 5: Visualize the Reachable Trajectory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   figure; hold on;
   for t = 1:length(R)
       for s = 1:length(R{t})
           % Project to (x1, x2) plane
           proj = R{t}(s).affineMap([1 0 0 0 0; 0 1 0 0 0], [0; 0]);
           proj.plot('b', 0.1);   % Blue with transparency
       end
   end
   xlabel('x_1'); ylabel('x_2');
   title('Neural ODE Reachable States (x_1 vs x_2)');

For a fixed-point attractor, you should see the initial uncertainty ball
**shrink** as the system converges to the attractor.

Reachability Methods for ODEs
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Method
     - Use Case
     - Notes
   * - ``direct``
     - Linear ODEs (dx/dt = Ax + Bu)
     - Uses matrix exponential. Exact and fast for linear systems.
   * - ``ode45``
     - Nonlinear ODEs
     - Numerical integration via CORA. Sound over-approximation.
   * - ``krylov``
     - High-dimensional linear ODEs
     - Krylov subspace approximation. Efficient for large state spaces.

For **linear** Neural ODEs, use ``direct`` for the best precision.
For **nonlinear** dynamics (sigmoid, tanh activations in the ODE), use ``ode45``
with CORA's Taylor/zonotope methods.

Source Files
^^^^^^^^^^^^

- `NN/NeuralODEs/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/NeuralODEs>`_
