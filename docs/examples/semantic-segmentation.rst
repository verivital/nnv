Semantic Segmentation
=====================

.. rst-class:: lead

   Verify pixel-level predictions of segmentation networks. This tutorial
   covers ImageStar construction for segmentation, handling large output
   dimensions, and scalability strategies.

----

What You Will Learn
-------------------

- Why semantic segmentation is harder to verify than classification
- How to create ImageStar input sets for segmentation networks
- How to use relaxed reachability for tractable analysis
- When to switch to probabilistic (conformal prediction) verification

Background
----------

Semantic segmentation networks assign a class label to **every pixel** in an
image. Unlike classification (which has a single output label), a segmentation
network for a 128x128 image with 11 classes produces 128 * 128 * 11 = 180,224
output values. This makes full reachability analysis extremely expensive.

NNV addresses this through:

1. **Relaxed reachability** (``relaxFactor > 0``) -- trade precision for speed
2. **Conformal prediction** -- probabilistic guarantees when exact methods are intractable

M2NIST Segmentation
---------------------

The M2NIST dataset places multiple MNIST digits on a canvas. The segmentation
network labels each pixel by digit identity.

.. code-block:: matlab

   % Load the segmentation network
   load('m2nist_seg_net.mat');
   net = matlab2nnv(seg_net);

   % Load a test image
   load('m2nist_test.mat');
   img = test_images(:,:,:,1);   % e.g., 64x84 single-channel image

   % Create L-inf perturbation set
   epsilon = 0.005;   % Small perturbation (segmentation is sensitive)
   lb = max(single(img) - epsilon, 0);
   ub = min(single(img) + epsilon, 1);
   IS = ImageStar(lb, ub);

   % Use relaxed reachability for tractability
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';
   reachOptions.relaxFactor = 0.5;   % Partial relaxation (0=precise, 1=fast)

   % Compute reachable output set
   tic;
   output_sets = net.reach(IS, reachOptions);
   fprintf('Segmentation reachability: %.1f seconds\n', toc);

   % Check per-pixel robustness
   % For each pixel, verify that the correct class has the highest score

Scalability
^^^^^^^^^^^

For large segmentation networks (e.g., full CamVid street scene models),
even relaxed reachability may be intractable. In these cases, use
:doc:`/user-guide/conformal-prediction` for probabilistic guarantees:

.. code-block:: matlab

   reachOptions.coverage = 0.99;
   reachOptions.confidence = 0.99;
   reachOptions.train_mode = 'Linear';
   result = verify_robustness_cp(net, IS, reachOptions, target, numClasses);

CamVid Street Scene Segmentation
----------------------------------

The CamVid dataset contains urban driving scenes with 11 classes (road,
building, car, pedestrian, etc.). NNV includes examples using dilated
convolution architectures with batch normalization.

These models demonstrate:

- Multi-scale feature learning with dilation factors (1, 2, 4)
- Batch normalization in the verification pipeline
- Handling very large output dimensions (per-pixel, per-class)

Source Files
^^^^^^^^^^^^

- `NN/SemanticSegmentation/M2NIST/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/SemanticSegmentation/M2NIST>`_
- `NN/SemanticSegmentation/Camvid/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/SemanticSegmentation>`_
