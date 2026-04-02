Custom Verification Workflows
=============================

.. rst-class:: lead

   End-to-end patterns for verifying your own neural network, from loading
   to result interpretation. Covers FFNNs, CNNs, GNNs, and NNCS.

----

Generic FFNN Workflow
---------------------

.. code-block:: matlab

   %% 1. Load or build network
   % From ONNX:
   dlnet = importNetworkFromONNX('model.onnx', InputDataFormats='BC');
   net = matlab2nnv(dlnet);
   % Or build manually:
   net = NN({FullyConnectedLayer(W1,b1), ReluLayer(), FullyConnectedLayer(W2,b2)});

   %% 2. Define input set
   input_set = Star(lb, ub);  % Box region

   %% 3. Configure and run reachability
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';
   output_sets = net.reach(input_set, reachOptions);

   %% 4. Check safety property
   unsafe_region = HalfSpace(G, g);  % Gx <= g defines unsafe outputs
   result = verify_specification(output_sets, unsafe_region);
   % 1=safe, 0=unsafe, 2=unknown

   %% 5. Extract output bounds
   [lb_out, ub_out] = output_sets(1).getRanges();

CNN with ImageStar
------------------

.. code-block:: matlab

   net = matlab2nnv(cnn_model);

   % L-inf perturbation around an image
   epsilon = 0.01;
   lb = max(single(img) - epsilon, 0);
   ub = min(single(img) + epsilon, 1);
   IS = ImageStar(lb, ub);

   reachOptions.reachMethod = 'approx-star';
   result = net.verify_robustness(IS, reachOptions, target_class);

GNN with GraphStar
------------------

.. code-block:: matlab

   gnn = GNN(layers, A_norm, adj_list, E);

   % Selective node feature perturbation
   eps_matrix = zeros(numNodes, numFeatures);
   eps_matrix(:, perturb_features) = epsilon;
   GS = GraphStar(X, -eps_matrix, eps_matrix);

   reachOptions.reachMethod = 'approx-star';
   output_sets = gnn.reach(GS, reachOptions);

   % Per-node safety check
   [lb, ub] = output_sets.getRanges();

NNCS Closed-Loop
----------------

.. code-block:: matlab

   controller = matlab2nnv(nn_controller);
   plant = NonLinearODE(6, 1, @dynamics, reachStep, controlPeriod, C);
   nncs = NonlinearNNCS(controller, plant);
   nncs.feedbackMap = [0];

   reachPRM.init_set = Star(x0_lb, x0_ub);
   reachPRM.ref_input = Star(ref_lb, ref_ub);
   reachPRM.numSteps = 50;
   reachPRM.reachMethod = 'approx-star';
   [R, rT] = nncs.reach(reachPRM);

Writing Custom Safety Specifications
--------------------------------------

Safety properties are defined as HalfSpace constraints on the output:

.. code-block:: matlab

   % Property: output y1 should be >= 0  (unsafe if y1 < 0)
   % Unsafe region: y1 <= 0  →  [-1 0 ... 0] * y <= 0
   unsafe = HalfSpace([-1 zeros(1, output_dim-1)], 0);

   % Property: class c should have highest score (robustness)
   % Unsafe: any other class j has score >= class c
   % For each j != c:  y_j - y_c >= 0  →  e_j - e_c (as row) * y <= 0
   for j = 1:num_classes
       if j ~= target
           G = zeros(1, num_classes);
           G(j) = 1; G(target) = -1;
           unsafe_j = HalfSpace(G, 0);
           result_j = verify_specification(output_sets, unsafe_j);
       end
   end

   % Combined: use verify_robustness() which handles this automatically
   result = net.verify_robustness(input_set, reachOptions, target);

Choosing the Right Method
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Scenario
     - Method
     - Why
   * - Small FFNN, need proof
     - ``exact-star``
     - Complete: can prove both safety and unsafety
   * - Medium CNN
     - ``approx-star``
     - Good balance of precision and speed (default)
   * - Large network, quick screen
     - ``approx-zono``
     - Fastest, no LP solving required
   * - Very large / intractable
     - Conformal prediction
     - Probabilistic guarantees, architecture-agnostic
   * - Neural ODE (linear)
     - ``direct``
     - Matrix exponential, exact for linear systems
   * - Neural ODE (nonlinear)
     - ``ode45``
     - CORA integration for nonlinear dynamics
