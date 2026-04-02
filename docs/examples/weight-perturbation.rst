Weight Perturbation (ModelStar)
================================

.. rst-class:: lead

   Verify that neural networks remain robust under weight perturbations
   caused by quantization, hardware faults, or model compression.

----

What You Will Learn
-------------------

- How to specify weight perturbations using ModelStar's perturbation maps
- How to verify a single layer's robustness against weight compression
- How to interpret results in the context of quantization bit-widths
- How ModelStar compares against existing weight perturbation bounds

Background
----------

Most verification research assumes fixed network weights and analyzes robustness
to input perturbations. But deployed networks face weight uncertainty from
quantization, compression, and hardware limitations. ModelStar addresses this
by incorporating weight perturbations directly into the Star set reachability
framework.

See :doc:`/theory/weight-perturbation` for the theoretical foundations.

Prerequisites
^^^^^^^^^^^^^

- NNV installed with the ``WPutils`` module
- A trained MATLAB neural network (or ONNX model converted via ``matlab2nnv``)

Step 1: Load the Network
--------------------------

.. code-block:: matlab

   % Load the MNIST MLP (5 hidden FC layers: 1024, 512, 256, 256, 256)
   load('mnist_mlp.mat');
   net = matlab2nnv(trained_net);

   % Select a test image
   load('mnist_test.mat');
   img = test_images(:,:,:,1);
   target = test_labels(1);

Step 2: Define Weight Perturbation
------------------------------------

Perturbation magnitudes are specified relative to each layer's **weight range**
(the difference between the maximum and minimum scalar weights in the layer).
This is motivated by the direct connection to hardware: the weight range
translates to the minimum tolerance threshold for programming weights to
memristors in CiM devices.

.. code-block:: matlab

   % Choose which layer to perturb
   layer_idx = 6;   % fc_6 (output layer, 2560 weights, 10 neurons)

   % Perturbation magnitude as fraction of weight range
   epsilon = 0.01;  % 1% of weight range

   % Apply perturbation to the entire weight matrix of this layer
   net = perturb_whole_layer(net, layer_idx, epsilon);

   % For selective perturbation (e.g., specific weight subsets):
   % net = add_pert(net, layer_idx, epsilon, perturbation_map);

Step 3: Verify Robustness
--------------------------

.. code-block:: matlab

   % Create singleton input (single image, no input perturbation)
   input_set = Star(single(img(:)), single(img(:)));

   % Verify
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';
   result = net.verify_robustness(input_set, reachOptions, target);

   switch result
       case 1, fprintf('ROBUST: Classification unchanged under weight perturbation.\n');
       case 0, fprintf('NOT ROBUST: Weight perturbation can change classification.\n');
       case 2, fprintf('UNKNOWN: Over-approximation too coarse to determine.\n');
   end

Connection to Quantization
---------------------------

The perturbation magnitudes used in ModelStar directly correspond to quantization
error levels. For an n-bit quantization scheme, the maximum relative error induced
by rounding (half of the least significant bit) is:

.. list-table::
   :header-rows: 1
   :widths: 20 30 30

   * - Bit Width
     - Max Relative Error
     - Approx. Perturbation Level
   * - 10-bit
     - 0.0488%
     - ~0.05%
   * - 9-bit
     - 0.0977%
     - ~0.1%
   * - 8-bit
     - 0.1953%
     - ~0.2%
   * - 7-bit
     - 0.3906%
     - ~0.4%

This means ModelStar's experimental perturbation levels (0.05%, 0.1%, 0.2%, 0.4%)
directly capture robustness against compression-induced quantization errors at
corresponding bit-widths.

Experimental Results
---------------------

ModelStar was evaluated on an MNIST MLP against two baselines:

- **Certificated-Robust** (Weng et al., AAAI 2020): linear bound propagation
- **Formal-Robust** (Tsai et al., NeurIPS 2021): pair-wise class margin bounds

.. raw:: html

   <div class="table-scroll">

.. list-table::
   :header-rows: 1

   * - Layer
     - Perturbation
     - ModelStar
     - Cert.-Robust
     - Formal-Robust
   * - fc_4
     - 0.1%
     - 100/100
     - 99/100
     - 0/100
   * - fc_4
     - 0.2%
     - 69/100
     - 13/100
     - 0/100
   * - fc_5
     - 0.1%
     - 100/100
     - 99/100
     - 89/100
   * - fc_5
     - 0.2%
     - 92/100
     - 41/100
     - 12/100
   * - fc_6
     - 0.5%
     - 100/100
     - 100/100
     - 100/100
   * - fc_6
     - 1.0%
     - 83/100
     - 83/100
     - 83/100

.. raw:: html

   </div>

At 0.2% perturbation on layer fc_4, ModelStar verifies **69/100 images** vs.
13 for Certificated-Robust -- a **56% improvement**.

Quantization Interpretation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From a quantization perspective, the results indicate that:

- All 100 images remain correctly classified under **7-bit weight compression**
  of layer fc_6 (max error ~0.39%, within the 0.5% level where 100% verified)
- All 100 images remain correctly classified under **9-bit compression** of
  layers fc_4 or fc_5 (max error ~0.1%)

Note that the current setup verifies compression of **one layer at a time**.
Multi-layer simultaneous perturbation is supported but yields lower verification
rates due to accumulated over-approximation.

Adapting to Your Own Network
-----------------------------

1. **Train your network** in MATLAB or load from ONNX
2. **Identify critical layers**: Later layers with fewer parameters are cheaper to verify
3. **Choose perturbation level**: Map your deployment's quantization bit-width to
   the corresponding error percentage
4. **Verify per-layer**: Start with the output layer and work backward
5. **Interpret results**: "Verified robust" means the classification is guaranteed
   correct under that perturbation level

Source Files
^^^^^^^^^^^^

- `Tutorial/NN/Collins_RUL/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NN/Collins_RUL>`_ -- RUL prediction with weight perturbation
