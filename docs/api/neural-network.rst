NN Class
========

.. rst-class:: lead

   The unified neural network class for feedforward, convolutional, recurrent,
   segmentation, and ODE-based architectures.

----

Constructors
------------

.. code-block:: matlab

   nn = NN(Layers)
   nn = NN(Layers, Connections)
   nn = NN(Layers, Connections, inputSize, outputSize, name)

- ``Layers`` -- cell array of layer objects (FullyConnectedLayer, ReluLayer, etc.)
- ``Connections`` -- table defining DAG connectivity (optional; omit for sequential networks)
- ``inputSize``, ``outputSize`` -- dimension arrays
- ``name`` -- string identifier

Properties
----------

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Property
     - Type
     - Description
   * - ``Layers``
     - cell array
     - Ordered list of layer objects
   * - ``Connections``
     - table
     - DAG connectivity (source/destination layer names)
   * - ``Name``
     - string
     - Network name
   * - ``numLayers``
     - int
     - Number of layers
   * - ``numNeurons``
     - int
     - Total neuron count across all layers
   * - ``InputSize``
     - array
     - Input dimensions
   * - ``OutputSize``
     - array
     - Output dimensions
   * - ``reachSet``
     - cell array
     - Computed reachable set per layer (after calling reach)
   * - ``reachTime``
     - array
     - Computation time per layer (after calling reach)

Methods
-------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Method
     - Description
   * - ``y = nn.evaluate(x)``
     - Forward pass on a single concrete input vector/image/tensor
   * - ``outputSet = nn.reach(inputSet, reachOptions)``
     - Compute reachable output set from an input set (Star, ImageStar, etc.)
   * - ``result = nn.verify_robustness(inputSet, reachOptions, target)``
     - Verify classification robustness. Returns: 1=robust, 0=not robust, 2=unknown.
       For exact methods, unknown is treated as not robust (returns 0).
   * - ``[result, X] = nn.verify_vnnlib(vnnlib_file, reachOptions)``
     - Verify properties from a VNNLIB file. Returns result (0/1/2) and the input set X.
   * - ``result = nn.verify_safety(inputSet, unsafeRegion, reachOptions)``
     - Verify safety against a HalfSpace unsafe region
   * - ``result = nn.verify_segmentation(inputSet, reachOptions, target)``
     - Verify pixel-level robustness for semantic segmentation networks
   * - ``label = nn.classify(inputValue, reachOptions)``
     - Classify input(s); handles both single points and sets
   * - ``rb = nn.checkRobust(outputSet, target)``
     - Check if a precomputed reachable set satisfies robustness for target class
   * - ``counter = nn.falsify(I, U, n_samples)``
     - Find counterexample inputs violating specification via random sampling
   * - ``result = nn.verify_sequence_robustness(input_seq, eps, target, reachOptions)``
     - Verify robustness for RNN/sequence classification over time steps

Example
-------

.. code-block:: matlab

   % Build a simple network
   layers = {FullyConnectedLayer(W1, b1), ReluLayer(), FullyConnectedLayer(W2, b2)};
   net = NN(layers);

   % Evaluate
   y = net.evaluate([0.5; 0.5]);

   % Verify
   input_set = Star([-1; -1], [1; 1]);
   reachOptions.reachMethod = 'approx-star';
   output_sets = net.reach(input_set, reachOptions);
   result = net.verify_robustness(input_set, reachOptions, 1);
