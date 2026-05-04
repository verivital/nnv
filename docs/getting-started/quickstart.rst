Quick Start
===========

.. rst-class:: lead

   Verify your installation, then run your first neural network
   verification in under five minutes.

----

Test Your Installation
----------------------

After installing NNV (see :doc:`installation`), verify everything works.

.. note::

   Run ``install.m`` on first use. If you restart MATLAB, run ``startup_nnv.m``
   to restore paths (or run ``savepath`` after installation to make it permanent).

.. code-block:: matlab

   cd examples/QuickStart
   test_installation      % Checks NNV version, Star class, layers, and NN creation
   simple_verification    % Runs a complete verification workflow

The ``test_installation`` script checks that:

- ``NNVVERSION()`` returns the expected version
- Star sets can be created from bounds
- FullyConnectedLayer and ReLU layers work correctly
- A simple neural network can be assembled and evaluated

Your First Verification
-----------------------

Here is the core NNV workflow -- define a network, create an input set,
compute reachable outputs, and check a safety property:

.. code-block:: matlab

   %% Step 1: Create a simple neural network (2 inputs → 4 hidden [ReLU] → 2 outputs)
   W1 = [1 0; 0 1; 1 1; -1 0];  % 4x2 weight matrix
   b1 = [0; 0; 0; 0];
   W2 = [1 0 1 0; 0 1 0 1];     % 2x4 weight matrix
   b2 = [0; 0];

   layers = {FullyConnectedLayer(W1, b1), ReLULayer(), FullyConnectedLayer(W2, b2)};
   net = NN(layers);

   %% Step 2: Define an input set (box region)
   lb = [-1; -1];   % lower bounds
   ub = [1; 1];     % upper bounds
   input_set = Star(lb, ub);

   %% Step 3: Compute reachable output set
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';
   output_sets = net.reach(input_set, reachOptions);

   %% Step 4: Analyze output bounds
   for i = 1:length(output_sets)
       [lb_out, ub_out] = output_sets(i).getRanges();
       fprintf('Output set %d bounds:\n', i);
       fprintf('  y1: [%.4f, %.4f]\n', lb_out(1), ub_out(1));
       fprintf('  y2: [%.4f, %.4f]\n', lb_out(2), ub_out(2));
   end

   %% Step 5: Check safety property (y1 >= -1 for all inputs)
   unsafe_region = HalfSpace([1 0], -1);  % y1 <= -1
   result = verify_specification(output_sets, unsafe_region);
   % result: 1 = safe (property holds), 0 = unsafe (violated), 2 = unknown

Key Concepts
------------

**Input sets**: Define the region of inputs you want to verify. Use ``Star(lb, ub)``
for box regions or construct more complex polytopic sets with constraints.

**Reachability**: ``net.reach()`` propagates the input set through every layer,
computing the set of all possible outputs. The method (``exact-star``,
``approx-star``, ``approx-zono``) controls the precision/speed tradeoff.

**Verification**: ``verify_specification()`` checks whether any output in the
reachable set violates a safety property defined as a HalfSpace constraint.

See :doc:`/user-guide/verification-methods` for a full guide to all methods.

Loading an ONNX Model
---------------------

Most real-world workflows start by loading a pre-trained ONNX model:

.. code-block:: matlab

   % Load ONNX network and convert to NNV format
   dlnet = importNetworkFromONNX('model.onnx', InputDataFormats='BC');
   net = matlab2nnv(dlnet);

   % Create input set from VNNLIB specification
   [lb, ub, prop] = load_vnnlib('property.vnnlib');
   input_set = Star(lb, ub);

   % Verify
   reachOptions.reachMethod = 'approx-star';
   res = net.verify_vnnlib('property.vnnlib', reachOptions);

See :doc:`/user-guide/onnx-vnnlib` for complete ONNX and VNNLIB workflow details.

Troubleshooting
---------------

If you encounter issues, run:

.. code-block:: matlab

   check_nnv_setup()

This diagnoses common problems including missing toolboxes, submodule issues,
and Python environment configuration.

Next Steps
----------

- :doc:`/user-guide/architectures` -- All supported network types and layers
- :doc:`/user-guide/set-representations` -- Understanding Star sets, ImageStars, and more
- :doc:`/user-guide/verification-methods` -- Choosing the right verification method
- :doc:`/examples/index` -- Worked examples across 10+ domains
