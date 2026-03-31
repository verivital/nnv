Control Systems
===============

.. rst-class:: lead

   Verify the safety of closed-loop control systems with neural network
   controllers. This tutorial covers adaptive cruise control, emergency
   braking, and inverted pendulum stabilization.

----

What You Will Learn
-------------------

- How to define plant dynamics (linear, nonlinear, discrete, continuous)
- How to compose a neural network controller with a plant model
- How to set initial conditions and run bounded-time reachability
- How to check safety specifications over time (e.g., minimum following distance)
- How to visualize reachable state trajectories

Adaptive Cruise Control (ACC)
------------------------------

The ACC system maintains a safe following distance between an ego vehicle
and a lead vehicle using a neural network controller. This is a 6-state
nonlinear system verified over 5 seconds of driving.

Prerequisites
^^^^^^^^^^^^^

- ``controller_5_20.mat`` -- trained NN controller (in the Tutorial folder)
- ``dynamicsACC.m`` -- plant dynamics function (in the same folder)

Step 1: Load the Neural Network Controller
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Load the trained controller (5 hidden layers, 20 neurons each)
   load('controller_5_20.mat');
   net = matlab2nnv(controller);

   % The controller maps 3 inputs to 1 output:
   %   Inputs:  [relative_distance, relative_velocity, ego_velocity]
   %   Output:  [acceleration_command]

Step 2: Define the Plant Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ACC plant is a 6-state nonlinear system modeling both vehicles:

.. code-block:: matlab

   % States: [x_lead, v_lead, a_lead, x_ego, v_ego, a_ego]
   %   x = position, v = velocity, a = acceleration
   %
   % The dynamics include friction effects and bounded acceleration

   reachStep = 0.01;       % ODE integration step (seconds)
   controlPeriod = 0.1;    % Controller sampling period (seconds)

   % Output matrix: maps 6 states to 3 controller inputs
   %   output = [x_lead - x_ego; v_lead - v_ego; v_ego]
   output_mat = [-1 0 0 1 0 0;    % relative distance
                  0 -1 0 0 1 0;    % relative velocity
                  0  0 0 0 1 0];   % ego velocity

   % Create nonlinear plant (wraps CORA for reachability)
   plant = NonLinearODE(6, 1, @dynamicsACC, reachStep, controlPeriod, output_mat);
   plant.options.taylorTerms = 4;       % Taylor expansion order
   plant.options.zonotopeOrder = 20;    % Zonotope complexity limit

.. admonition:: What is NonLinearODE?

   ``NonLinearODE(nStates, nInputs, @dynamics, reachStep, controlPeriod, C)``
   creates a continuous-time nonlinear plant model. NNV uses the CORA toolbox
   internally to compute reachable sets of the ODE using Taylor expansion and
   zonotope propagation. The ``reachStep`` controls the ODE integration
   granularity, while ``controlPeriod`` is the interval between controller
   actions. See :doc:`/user-guide/nncs` for all plant model types.

Step 3: Compose the NNCS
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Create the closed-loop system
   nncs = NonlinearNNCS(net, plant);

   % feedbackMap defines which plant outputs feed into the controller
   % [0] means the controller receives the plant output directly
   nncs.feedbackMap = [0];

Step 4: Define Initial Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Initial state ranges for both vehicles
   %          [x_lead,  v_lead,  a_lead,  x_ego,  v_ego,  a_ego]
   lb_state = [  90;      29;       0;      10;     30;      0  ];
   ub_state = [ 110;      30;       0;      11;     30.5;    0  ];
   init_set = Star(lb_state, ub_state);

   % Reference inputs to the controller
   target_speed = 30;    % m/s
   time_gap = 1.4;       % seconds (desired following time gap)
   ref_input = Star([target_speed; time_gap], [target_speed; time_gap]);

Step 5: Run Reachability Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Configure reachability parameters
   reachPRM = struct;
   reachPRM.init_set = init_set;
   reachPRM.ref_input = ref_input;
   reachPRM.numSteps = 50;               % 50 control steps = 5 seconds
   reachPRM.reachMethod = 'approx-star';
   reachPRM.numCores = 1;                % Set > 1 for parallel computation

   % Compute reachable states over time
   tic;
   [R, rT] = nncs.reach(reachPRM);
   fprintf('Reachability completed in %.1f seconds\n', toc);

   % R contains the reachable state sets at each time step
   % R{t} is an array of Star sets representing all possible states at step t

Step 6: Check Safety
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Safety specification: relative distance > time_gap * ego_velocity + buffer
   % Extract relative distance at each time step
   safe = true;
   for t = 1:length(R)
       for s = 1:length(R{t})
           % Project to relative distance (x_lead - x_ego)
           dist_map = [-1 0 0 1 0 0];   % x_lead - x_ego
           dist_set = R{t}(s).affineMap(dist_map, 0);
           [d_lb, d_ub] = dist_set.getRanges();

           % Check minimum safe distance (simplified: > 10 meters)
           if d_lb < 10
               fprintf('Potential safety violation at step %d: min distance = %.2f m\n', t, d_lb);
               safe = false;
           end
       end
   end

   if safe
       fprintf('SAFE: Minimum following distance maintained for all 5 seconds.\n');
   end

Automated Emergency Braking System (AEBS)
------------------------------------------

The AEBS uses a reinforcement learning controller with saturation activation
functions. The 3-state system (distance, velocity, acceleration) uses a discrete
linear plant.

.. code-block:: matlab

   % Load RL controller and transformer networks
   load('models/controller.mat');   % RL controller weights
   load('models/transform.mat');    % Output transformer weights

   % Build the controller as a chain of NNV layers:
   %   1. Normalizer (scales inputs to [-1, 1])
   %   2. RL controller (3 layers with satlin activation)
   %   3. Transformer (maps controller output to brake command)
   %   4. Output scaler (scales to [-500, 500] N)

   % Define 3-state discrete linear plant
   % State: [distance (m), velocity (m/s), acceleration (m/s^2)]
   plant = DLinearODE(A, B, C, D);   % A, B, C, D are discrete-time matrices

   % Compose and verify over 50 control steps
   nncs = DLinearNNCS(controller, plant);
   nncs.feedbackMap = [0];

The AEBS example demonstrates how to handle multi-component controller
architectures where normalizer, controller, and output transformer are
separate neural network blocks composed together.

**Source:** `Tutorial/NNCS/AEBS/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NNCS/AEBS>`_

Inverted Pendulum
------------------

A classic control benchmark: a neural network stabilizes an inverted pendulum
using state feedback from a 4-state discrete-time linear plant.

.. code-block:: matlab

   % Load controller and plant
   load('controller.mat');   % NN weights: W1, b1, W2, b2
   load('sysd.mat');         % Discrete system: A, B, C, D, Ts

   % Build 2-layer NN controller
   L1 = FullyConnectedLayer(W1, b1);
   L2 = ReluLayer();         % poslin activation
   L3 = FullyConnectedLayer(W2, b2);
   controller = NN({L1, L2, L3});

   % Create discrete linear plant and compose NNCS
   plant = DLinearODE(A, B, C, D);
   ncs = DLinearNNCS(controller, plant);
   ncs.feedbackMap = [0];   % Full state feedback

   % Initial conditions
   %   x1: position [-0.1, -0.05] m
   %   x2: angle [0.85, 0.9] rad (near vertical)
   %   x3, x4: velocities = 0
   lb = [-0.1; 0.85; 0; 0];
   ub = [-0.05; 0.9; 0; 0];
   init_set = Star(lb, ub);

   % Reachability over 20 steps (2 seconds)
   reachPRM.init_set = init_set;
   reachPRM.numSteps = 20;
   reachPRM.reachMethod = 'approx-star';
   [R, rT] = ncs.reach(reachPRM);

   % Visualize position vs angle trajectory
   maps = [1 0 0 0; 0 1 0 0];   % Extract position and angle
   figure; hold on;
   for t = 1:length(R)
       for s = 1:length(R{t})
           proj = R{t}(s).affineMap(maps, [0; 0]);
           proj.plot();           % Plot 2D projection
       end
   end
   xlabel('Position (m)');
   ylabel('Angle (rad)');
   title('Inverted Pendulum Reachable States');

**Source:** `Tutorial/NNCS/InvertedPendulum/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NNCS/InvertedPendulum>`_

Key NNCS Parameters
--------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``controlPeriod``
     - Time between controller actions (seconds). The plant evolves continuously
       between control inputs.
   * - ``reachStep``
     - ODE integration step within each control period. Smaller values give
       more precise plant reachability but cost more computation.
   * - ``numSteps``
     - Number of control periods to simulate. Total time = numSteps * controlPeriod.
   * - ``feedbackMap``
     - Defines which plant outputs are fed back to the controller. ``[0]`` means
       the controller receives the full plant output vector.
   * - ``taylorTerms``
     - Taylor expansion order for nonlinear plant reachability (CORA). Higher
       values are more precise but slower. Typical: 4-6.
   * - ``zonotopeOrder``
     - Maximum zonotope generator count for plant reachability. Higher values
       allow more precise representation. Typical: 20-50.

Source Files
^^^^^^^^^^^^

- `Tutorial/NNCS/ACC/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NNCS/ACC>`_
- `Tutorial/NNCS/AEBS/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NNCS/AEBS>`_
- `Tutorial/NNCS/InvertedPendulum/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NNCS/InvertedPendulum>`_
- `NN/dcdc_buck/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/dcdc_buck>`_
