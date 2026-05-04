Architecture Overview
=====================

.. rst-class:: lead

   How NNV is structured internally: the two-module pipeline, directory
   layout, and how networks compose layers for reachability.

----

Two-Module Pipeline
-------------------

NNV is composed of two primary modules:

1. **Computation Engine**: Parses models (ONNX, MATLAB), constructs NNV
   representations, and performs layer-by-layer reachability analysis using
   set-based abstractions (Star, ImageStar, VolumeStar, GraphStar).

2. **Analyzer**: Consumes the reachable sets to verify safety properties,
   generate counterexamples, and visualize results.

Directory Structure
-------------------

.. code-block:: text

   code/nnv/engine/
   ├── nn/                    # Neural network classes
   │   ├── NN.m               # Main unified NN class
   │   ├── GNN.m              # Graph neural network class
   │   ├── layers/            # 48 layer type implementations
   │   │   ├── FullyConnectedLayer.m    # Canonical affine layer (~650 lines)
   │   │   ├── ReluLayer.m              # ReLU with exact/approx reachability
   │   │   ├── Conv2DLayer.m            # 2D convolution
   │   │   ├── GCNLayer.m              # Graph convolution
   │   │   └── ...
   │   ├── funcs/             # Activation function math (PosLin, LogSig, etc.)
   │   └── Prob_reach/        # Probabilistic verification (Python trainers)
   ├── set/                   # Set representations
   │   ├── Star.m             # ~2050 lines, the fundamental set type
   │   ├── ImageStar.m, VolumeStar.m, GraphStar.m
   │   ├── Zono.m, ImageZono.m, Box.m, HalfSpace.m
   │   └── Conversion.m, SetTree.m
   ├── nncs/                  # Control system classes
   │   ├── LinearODE.m, NonLinearODE.m, HybridA.m
   │   └── LinearNNCS.m, NonlinearNNCS.m, ...
   ├── utils/                 # Utilities
   │   ├── matlab2nnv.m       # MATLAB network → NNV conversion
   │   ├── onnx2nnv.m, load_vnnlib.m
   │   ├── verify_specification.m, lpsolver.m
   │   └── WPutils.m, Reduction.m, ...
   ├── cora/                  # CORA toolbox (submodule)
   └── hyst/                  # HyST hybrid systems tool (submodule)

How NN Composes Layers
----------------------

The ``NN`` class stores layers in a cell array and optionally a ``Connections``
table for DAG (skip-connection) networks. The main ``reach()`` method dispatches
based on whether connections exist:

.. code-block:: matlab

   function outputSet = reach(obj, inputSet, reachOptions)
       % ... parse reachOptions into obj properties ...
       if isempty(obj.Connections)
           outputSet = obj.reach_noConns(inputSet);
       else
           outputSet = obj.reach_withConns(inputSet);
       end
   end

**Sequential networks** (``reach_noConns``):

The core loop iterates through layers, passing 6 arguments to each layer's
``reach()`` via ``varargin``:

.. code-block:: matlab

   function outSet = reach_noConns(obj, inSet)
       rs = inSet;
       for i = 1:obj.numLayers
           rs = obj.Layers{i}.reach(rs, obj.reachMethod, obj.reachOption, ...
                                     obj.relaxFactor, obj.dis_opt, obj.lp_solver);
           obj.reachSet{i} = rs;    % store for inspection
           obj.reachTime(i) = toc;  % record timing
       end
       outSet = rs;
   end

The 6 arguments passed to every ``layer.reach()`` are:

1. ``inputSet`` -- Star, ImageStar, VolumeStar, or array of Star sets
2. ``method`` -- string: ``'approx-star'``, ``'exact-star'``, ``'abs-dom'``, ``'relax-star-0.5'``, ``'approx-zono'``
3. ``option`` -- string: ``'single'`` or ``'parallel'`` (set when ``numCores > 1``)
4. ``relaxFactor`` -- numeric: controls LP vs interval bound fraction for approx-star
5. ``dis_opt`` -- display option: ``'display'`` or ``[]``
6. ``lp_solver`` -- string: ``'linprog'`` or ``'glpk'``

.. note::

   Not all layers use all 6 arguments. Affine layers (FullyConnectedLayer,
   Conv2DLayer) only use the input set, method, and option -- they ignore
   relaxFactor, dis_opt, and lp_solver since their reachability is exact.
   Nonlinear layers (ReluLayer) use all 6.

**DAG networks** (``reach_withConns``):

Uses the ``Connections`` table to route sets between named layers. Handles
multi-input destinations (e.g., AdditionLayer receiving from two branches)
by parsing destination specifications like ``'add_1/in1'``, ``'add_1/in2'``.

How GNN Differs from NN
------------------------

The ``GNN`` class overrides the reach loop because graph layers require
additional arguments (adjacency matrix, edge features) that standard layers
don't accept:

.. code-block:: matlab

   % GNN.reach() - simplified from actual code
   for i = 1:obj.numLayers
       layer = obj.Layers{i};
       if isa(layer, 'GCNLayer')
           % GCN needs normalized adjacency matrix as 3rd argument
           rs = layer.reach(rs, obj.A_norm, obj.reachMethod, obj.reachOption);
       elseif isa(layer, 'GINELayer')
           % GINE needs edge features, edge list, and edge weights
           rs = layer.reach(rs, obj.E, obj.adj_list, obj.reachMethod, ...
                            obj.reachOption, obj.relaxFactor, obj.lp_solver, ...
                            obj.edge_weights);
       elseif isa(layer, 'ReluLayer')
           % Activation layers: convert GraphStar → Star → apply → convert back
           numNodes = rs.numNodes;
           numFeatures = rs.numFeatures;
           S = rs.toStar();                  % Flatten GraphStar to Star
           S_out = layer.reach(S, obj.reachMethod, obj.reachOption, ...
                               obj.relaxFactor, obj.dis_opt, obj.lp_solver);
           % Reshape Star back to GraphStar
           new_V = reshape(S_out.V, [numNodes, numFeatures, S_out.nVar + 1]);
           rs = GraphStar(new_V, S_out.C, S_out.d, ...
                          S_out.predicate_lb, S_out.predicate_ub);
       else
           % Standard layers (e.g., FullyConnectedLayer after global pooling)
           rs = layer.reach(rs, obj.reachMethod, obj.reachOption, ...
                            obj.relaxFactor, obj.dis_opt, obj.lp_solver);
       end
   end

The key pattern is the **GraphStar ↔ Star conversion** at activation layers:
ReLU and other activations operate element-wise on flattened vectors, so the
graph structure is temporarily removed via ``toStar()`` and restored after
the activation by reshaping the Star's basis matrix ``V`` back to
``(numNodes, numFeatures, numBasis)`` format.

How NNCS Composes Controller + Plant
-------------------------------------

A neural network control system alternates between:

1. **Controller step**: Compute NN reachable output (control input) from plant feedback
2. **Plant step**: Propagate the control input through plant dynamics for one period
3. **Feedback**: Map plant output back to controller input via ``feedbackMap``

.. code-block:: matlab

   for t = 1:numSteps
       % Controller: NN reachability
       u_set = controller.reach(feedback_set, reachOptions);

       % Plant: ODE reachability (LinearODE uses matrix exponential,
       %   NonLinearODE uses CORA's Taylor/zonotope methods)
       state_set = plant.stepReach(u_set, ...);

       % Feedback mapping
       feedback_set = state_set.affineMap(feedbackMap, ...);
   end

The ``feedbackMap`` property defines which plant outputs are fed back to the
controller. ``[0]`` means the controller receives the plant output directly.

How Set Types Flow Through the Pipeline
-----------------------------------------

The set type determines what kind of network can be verified:

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Input Set Type
     - Network Types
     - Set Flow
   * - Star
     - FFNN
     - Star → FC (affineMap) → Star → ReLU (split/relax) → Star(s) → ...
   * - ImageStar
     - CNN
     - ImageStar → Conv2D (conv on generators) → ImageStar → ReLU → ... → Flatten → Star → FC → Star
   * - VolumeStar
     - 3D CNN
     - Same as ImageStar but 4D tensors through Conv3D/Pool3D
   * - GraphStar
     - GNN
     - GraphStar → GCN (A_norm * V * W) → GraphStar → ReLU (→Star→ReLU→GraphStar) → ...
   * - Star (in NNCS)
     - Controller + Plant
     - Star → NN reach → Star → ODE reach → Star (per time step)
