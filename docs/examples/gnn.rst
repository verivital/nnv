Graph Neural Networks
=====================

.. rst-class:: lead

   Verify GNN predictions for power grid safety under input uncertainty.
   This tutorial walks through GraphStar construction, GCN/GINE layer
   setup, and per-node voltage safety verification on IEEE bus systems.

----

What You Will Learn
-------------------

- How power systems are represented as graphs (buses=nodes, lines=edges)
- How to construct a GraphStar with selective feature perturbation
- How to create GCN and GINE layers from trained model weights
- How to verify per-node voltage safety specifications
- How to compare robustness across GNN architectures

Background
----------

Graph neural networks are used as fast surrogates for power flow solvers.
Given node features (active/reactive power injections) and edge features
(line impedances), a GNN predicts voltage magnitudes and angles at each bus.
Safety requires that predicted voltages stay within **[0.95, 1.05] per unit** --
violations could damage equipment or cause blackouts.

NNV verifies that small input perturbations (representing load uncertainty
or measurement noise) cannot push predicted voltages outside safe bounds.

Prerequisites
^^^^^^^^^^^^^

- Trained GCN and GINE models: ``models/gcn_ieee24.mat``, ``models/gine_ieee24.mat``
- Located in ``examples/NN/GNN/CAV26/``

Step 1: Load Model and Graph Data
-----------------------------------

.. code-block:: matlab

   % Load pre-trained GCN model for IEEE 24-bus system
   model_data = load('models/gcn_ieee24.mat');

   % Extract components:
   %   model_data.weights  -- layer weight matrices
   %   model_data.A_norm   -- normalized adjacency matrix (24x24)
   %   model_data.X_test   -- test node features (24 x num_features)
   %   model_data.Y_test   -- test labels (24 x num_outputs)
   %   model_data.mean_X, model_data.std_X -- normalization parameters

   A_norm = model_data.A_norm;    % Normalized adjacency matrix
   X = model_data.X_test(:,:,1);  % First test graph (24 nodes x F features)
   Y = model_data.Y_test(:,:,1);  % Ground truth voltages

   numNodes = size(X, 1);         % 24 buses
   fprintf('IEEE 24-bus system: %d nodes, %d features per node\n', ...
       numNodes, size(X, 2));

Step 2: Build the GNN in NNV
------------------------------

.. code-block:: matlab

   % Extract weight matrices for each GCN layer
   W1 = model_data.weights{1};   % First GCN layer weights
   b1 = model_data.biases{1};    % First GCN layer bias
   W2 = model_data.weights{2};
   b2 = model_data.biases{2};
   W3 = model_data.weights{3};
   b3 = model_data.biases{3};

   % Create NNV GCN layers
   L1 = GCNLayer(W1, b1);
   L2 = ReluLayer();
   L3 = GCNLayer(W2, b2);
   L4 = ReluLayer();
   L5 = GCNLayer(W3, b3);

   % Create GNN with graph structure
   layers = {L1, L2, L3, L4, L5};
   gnn = GNN(layers, A_norm);

Step 3: Create GraphStar with Selective Perturbation
------------------------------------------------------

Not all node features should be perturbed equally. For power systems, we
perturb only the **active and reactive power** features (columns 1 and 2),
not voltage or angle features:

.. code-block:: matlab

   epsilon = 0.005;   % Perturbation magnitude

   % Compute per-feature perturbation bounds
   % Only perturb features 1 (active power) and 2 (reactive power)
   perturb_features = [1, 2];
   range_per_col = max(X) - min(X);  % Feature range for scaling

   eps_matrix = zeros(numNodes, size(X, 2));
   for f = perturb_features
       eps_matrix(:, f) = range_per_col(f) * epsilon;
   end

   % Create GraphStar: nominal features ± scaled perturbation
   GS = GraphStar(X, -eps_matrix, eps_matrix);

.. admonition:: What is a GraphStar?

   A ``GraphStar(X, LB, UB)`` represents the set of all node feature matrices
   where each entry ``X(i,f)`` can vary between ``X(i,f) + LB(i,f)`` and
   ``X(i,f) + UB(i,f)``. It is the graph-structured analogue of ImageStar.
   See :doc:`/user-guide/set-representations` for details.

Step 4: Run Reachability
--------------------------

.. code-block:: matlab

   % Configure verification
   reachOpts = struct;
   reachOpts.reachMethod = 'approx-star';

   % Compute reachable output set
   tic;
   output_sets = gnn.reach(GS, reachOpts);
   fprintf('GNN reachability completed in %.2f seconds\n', toc);

   % output_sets contains a Star set for each node's predicted voltage

Step 5: Verify Per-Node Voltage Safety
-----------------------------------------

.. code-block:: matlab

   % Safety specification: voltage in [0.95, 1.05] per unit
   v_min = 0.95;
   v_max = 1.05;

   % Normalize specification to match model's output normalization
   v_min_norm = (v_min - model_data.mean_Y) / model_data.std_Y;
   v_max_norm = (v_max - model_data.mean_Y) / model_data.std_Y;

   % Check each bus
   results = zeros(numNodes, 1);
   for node = 1:numNodes
       % Extract this node's output set
       Y_node = output_sets(node);
       [lb, ub] = Y_node.getRanges();

       if lb >= v_min_norm && ub <= v_max_norm
           results(node) = 1;    % Verified safe
           status = 'SAFE';
       elseif ub < v_min_norm || lb > v_max_norm
           results(node) = 0;    % Violated (output outside safe bounds entirely)
           status = 'VIOLATED';
       else
           results(node) = 2;    % Unknown (bounds overlap with safety boundary)
           status = 'UNKNOWN';
       end

       fprintf('  Bus %2d: [%.4f, %.4f] -> %s\n', node, lb, ub, status);
   end

   fprintf('\nSummary: %d SAFE, %d VIOLATED, %d UNKNOWN (of %d buses)\n', ...
       sum(results==1), sum(results==0), sum(results==2), numNodes);

GCN vs GINE Comparison
------------------------

GINE layers incorporate edge features (line impedances), which can provide
more robust predictions. Run the same analysis with a GINE model:

.. code-block:: matlab

   % GINE layers use edge list and edge features instead of adjacency matrix
   model_gine = load('models/gine_ieee24.mat');

   % Create GINE layers
   L1 = GINELayer(W1, b1);   % Takes node AND edge features
   % ... (same pattern as GCN)

   % GNN with edge features
   gnn_gine = GNN(layers, model_gine.adj_list, model_gine.E);

   % For GINE+Edge perturbation: perturb edge features too
   eps_edge = 0.001;   % Fixed edge perturbation for impedance uncertainty
   GS_edge = GraphStar(X, -eps_matrix, eps_matrix, E, -eps_edge, eps_edge);

Typical results show GINE is more robust than GCN at larger perturbation
budgets, while GCN degrades faster. GINE with edge perturbation shows the
combined effect of node and edge uncertainty.

Interpreting Results
---------------------

For power systems, the verification tells operators:

- **SAFE buses**: Voltage predictions are guaranteed within safe limits
  despite load uncertainty -- no operator action needed
- **UNKNOWN buses**: The GNN might produce unsafe predictions under
  worst-case conditions -- these buses need monitoring
- **VIOLATED buses**: The GNN fails to guarantee safety under this
  uncertainty level -- consider retraining or using a more robust model

Experimental Results Summary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GNNV evaluation across IEEE systems and graph classification benchmarks:

- **High verification rates**: Near-perfect robustness for OPF across all
  systems; 70--99% for PF depending on perturbation level
- **Edge perturbation impact minimal**: Adding edge uncertainty (1% line
  parameter deviation) reduces robustness by at most 0.2%, with runtime
  overhead of 1.1--2.9x
- **Tighter than CORA**: GNNV consistently verifies more graphs than CORA's
  polynomial-zonotope abstractions -- up to 21.6% more at larger perturbation
  budgets (e.g., 99/100 vs 86/100 on IEEE-24 CFA at epsilon=0.005)
- **Scalable**: Subgraph verification completes in under a second for small
  perturbations, and under 90 seconds for IEEE-118 at epsilon=0.01

Helper Functions
^^^^^^^^^^^^^^^^^

The GNN examples provide convenience functions:

.. code-block:: matlab

   results = run_pf_verification('IEEE24', 'GCN', 0.01);
   results = run_pf_verification('IEEE39', 'GINE', 0.01);
   results = run_opf_verification('IEEE24', 'GINE', 0.01);

IEEE Bus Systems
^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 25 25 25

   * - System
     - Nodes
     - Edges
   * - IEEE 24-bus
     - 24
     - 92
   * - IEEE 39-bus
     - 39
     - 131
   * - IEEE 118-bus
     - 118
     - 476

Verification Status Codes
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Code
     - Meaning
   * - 1
     - Verified safe (output within specification bounds)
   * - 0
     - Violated (output outside specification bounds)
   * - 2
     - Unknown boundary (over-approximation overlaps with spec)
   * - 3
     - Unknown timeout (verification did not complete)
   * - -1
     - N/A (non-voltage bus, no specification applies)

Importing Your Own Models
^^^^^^^^^^^^^^^^^^^^^^^^^^

Models are trained in Python (PyTorch Geometric) and exported to MATLAB via
``.mat`` files. Note that **PyTorch uses (F_out x F_in) weight layout** --
transpose is required when loading into NNV's GCNLayer/GINELayer.

.. warning::

   GINELayer uses linear projections (not MLPs) for verification soundness.
   This means the NNV GINE implementation may differ slightly from the full
   PyTorch Geometric GINEConv with MLP update functions.

Source Files
^^^^^^^^^^^^

- `NN/GNN/CAV26/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/GNN/CAV26>`_
- `NN/GNN/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/GNN>`_
