Image Classification
====================

.. rst-class:: lead

   Verify that adversarial perturbations cannot fool image classifiers.
   This tutorial walks through complete MNIST and GTSRB verification
   workflows from start to finish.

----

What You Will Learn
-------------------

- How to load a pre-trained image classifier and convert it to NNV format
- How to define an adversarial perturbation region using ImageStar
- How to run robustness verification with ``verify_robustness()``
- How to interpret verification results (robust / not robust / unknown)
- How to visualize output ranges to understand classifier behavior

MNIST Digit Classification
----------------------------

This example verifies that small pixel perturbations to a handwritten digit
image cannot cause the classifier to change its prediction.

Prerequisites
^^^^^^^^^^^^^

- NNV installed and on the MATLAB path (``install.m`` has been run)
- Pre-trained MNIST model: ``mnist_model.mat`` (included in the Tutorial folder)
- MATLAB's built-in DigitDataset (ships with Deep Learning Toolbox, no download needed)

Step 1: Load and Convert the Network
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Load the pre-trained MNIST network
   mnist_model = load('mnist_model.mat');

   % Convert from MATLAB deep learning format to NNV format
   % matlab2nnv() handles SeriesNetwork, DAGNetwork, and dlnetwork
   net = matlab2nnv(mnist_model.net);

The ``matlab2nnv()`` function inspects each layer of the MATLAB network and
creates the corresponding NNV layer objects (FullyConnectedLayer, ReluLayer,
Conv2DLayer, etc.). After conversion, ``net`` is an NNV ``NN`` object that
supports reachability analysis.

Step 2: Load a Test Image
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Load MATLAB's built-in digit dataset (no download needed)
   digitDatasetPath = fullfile(toolboxdir('nnet'), ...
       'nndemos', 'nndatasets', 'DigitDataset');
   imds = imageDatastore(digitDatasetPath, ...
       'IncludeSubfolders', true, 'LabelSource', 'foldernames');

   % Read a specific test image
   idx = 10;                                    % Image index
   [img, fileInfo] = readimage(imds, idx);
   target = single(fileInfo.Label);             % True class label

   % img is a 28x28 grayscale image with pixel values in [0, 255]

Step 3: Create the Adversarial Input Set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An **ImageStar** represents the set of all images within a perturbation region.
For L-infinity robustness, we define a box around the original image where
each pixel can vary by a fixed amount:

.. code-block:: matlab

   % Define perturbation magnitude
   epsilon = 1;  % Each pixel can change by ±1 (out of [0, 255])

   % Compute perturbed bounds, clipped to valid pixel range
   lb = max(single(img) - epsilon, 0);       % Lower bound (no negative pixels)
   ub = min(single(img) + epsilon, 255);     % Upper bound (no overflow)

   % Create the ImageStar input set
   % ImageStar(lb, ub) represents ALL images where each pixel is in [lb, ub]
   IS = ImageStar(lb, ub);

.. admonition:: What is an ImageStar?

   An ImageStar is a set of images: every image where pixel (h,w,c) falls
   between ``lb(h,w,c)`` and ``ub(h,w,c)`` is contained in the set.
   NNV propagates this entire set through the network in one pass, computing
   the set of all possible outputs. See :doc:`/user-guide/set-representations`
   for the mathematical definition.

Step 4: Quick Falsification Check
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before running the expensive reachability analysis, check whether the network
gives the correct answer on the original image and the extreme corners of the
perturbation region:

.. code-block:: matlab

   % Evaluate on the original image
   y_orig = net.evaluate(single(img));
   [~, pred_orig] = max(y_orig);

   % Evaluate on the lower and upper bounds
   y_lb = net.evaluate(lb);
   [~, pred_lb] = max(y_lb);

   y_ub = net.evaluate(ub);
   [~, pred_ub] = max(y_ub);

   % If any prediction differs from the target, the network is NOT robust
   if pred_orig ~= target || pred_lb ~= target || pred_ub ~= target
       fprintf('Network is NOT robust: classification changes at bounds.\n');
   else
       fprintf('Bounds check passed. Running full verification...\n');
   end

This quick check catches many non-robust cases in milliseconds, saving the
cost of full reachability analysis.

Step 5: Run Robustness Verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Configure reachability options
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';  % Sound over-approximation

   % Run verification
   % Returns: 1 = robust, 0 = not robust, 2 = unknown
   result = net.verify_robustness(IS, reachOptions, target);

   % Interpret the result
   switch result
       case 1
           fprintf('VERIFIED: Network is robust for ALL images in the input set.\n');
       case 0
           fprintf('NOT ROBUST: A counterexample exists within the perturbation.\n');
       case 2
           fprintf('UNKNOWN: Over-approximation is too coarse to determine.\n');
           fprintf('Try: smaller epsilon, exact-star method, or relaxFactor=0.\n');
   end

.. admonition:: Understanding Verification Results

   - **1 (Robust)**: The correct class has the highest score for *every* image
     in the perturbation set. This is a formal guarantee.
   - **0 (Not Robust)**: There exists at least one image in the perturbation set
     that is classified differently. The network is vulnerable.
   - **2 (Unknown)**: The approximate method cannot determine robustness.
     The over-approximation of the output set overlaps between classes.
     This does NOT mean the network is vulnerable -- it means the analysis
     is inconclusive. Try a more precise method.

Step 6: Visualize Output Ranges
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   % Compute the reachable output set
   output_sets = net.reach(IS, reachOptions);

   % Extract output bounds for each class
   [lb_out, ub_out] = output_sets.getRanges();

   % Plot output ranges for all 10 digit classes
   figure;
   n_classes = length(lb_out);
   for i = 1:n_classes
       plot([i i], [lb_out(i) ub_out(i)], 'b-', 'LineWidth', 2); hold on;
       plot(i, y_orig(i), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
   end
   xlabel('Class (0-9)');
   ylabel('Network output score');
   title('Output ranges under adversarial perturbation');
   legend('Reachable range', 'Original output');

The blue bars show the range of possible output scores for each class across
all images in the perturbation set. If the bar for the correct class is entirely
above all other bars, the network is verified robust.

GTSRB Traffic Sign Recognition
--------------------------------

The same workflow applies to traffic sign classification, with minor differences
for image dimensions and class count.

.. code-block:: matlab

   % Load pre-trained GTSRB model
   gtsrb_model = load('gtsrb_model.mat');
   net = matlab2nnv(gtsrb_model.net);

   % Load dataset from NNV's data directory
   datapath = fullfile(nnvroot(), 'data', 'GTSRB');
   imds = imageDatastore(datapath, ...
       'IncludeSubfolders', true, 'LabelSource', 'foldernames');

   % Read and resize image (GTSRB model expects 30x29 pixels)
   [img, fileInfo] = readimage(imds, 1);
   img = imresize(img, [30 29]);
   target = single(fileInfo.Label);

   % Create input set and verify (same as MNIST)
   epsilon = 1;
   lb = max(single(img) - epsilon, 0);
   ub = min(single(img) + epsilon, 255);
   IS = ImageStar(lb, ub);

   reachOptions.reachMethod = 'approx-star';
   result = net.verify_robustness(IS, reachOptions, target);

.. note::

   The key difference from MNIST is image resizing: ``imresize(img, [30 29])``
   matches the network's expected input dimensions. The number of output classes
   is determined by ``net.OutputSize`` (43 for GTSRB).

MedMNIST and CIFAR-10
-----------------------

For larger networks (MedMNIST organ classifiers, CIFAR-10 CNNs), the same
workflow applies but may require:

- **Higher relaxFactor** (e.g., 0.5) for faster but coarser analysis
- **Conformal prediction** for networks where full reachability is intractable

See :doc:`medical-imaging` for MedMNIST examples and
:doc:`/user-guide/conformal-prediction` for the probabilistic verification setup.

Adapting to Your Own Classifier
--------------------------------

To verify your own image classifier:

1. **Train or load your network** in MATLAB (or import from ONNX)
2. **Convert to NNV**: ``net = matlab2nnv(your_network)``
3. **Choose epsilon**: Start small (e.g., 1/255) and increase
4. **Create ImageStar**: ``IS = ImageStar(max(img-eps, 0), min(img+eps, 1))``
   (normalize bounds to match your network's expected input range)
5. **Verify**: ``result = net.verify_robustness(IS, reachOptions, target)``

.. tip::

   If verification returns UNKNOWN for most samples, try:

   - Reducing epsilon (smaller perturbation region)
   - Using ``reachOptions.relaxFactor = 0`` (full LP optimization)
   - Using ``reachOptions.reachMethod = 'exact-star'`` (complete but slower)
   - Using ``reachOptions.numCores = 4`` (parallel computation)

Source Files
^^^^^^^^^^^^

- `Tutorial/NN/MNIST/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NN/MNIST>`_
- `Tutorial/NN/GTSRB/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NN/GTSRB>`_
- `NN/medmnist/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/medmnist>`_
- `NN/cifar10/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/cifar10>`_
