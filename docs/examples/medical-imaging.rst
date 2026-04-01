Medical Imaging
===============

.. rst-class:: lead

   Verify robustness of medical image classifiers and segmentation networks
   against clinically relevant perturbations. Based on the SPIE 2025
   tutorial materials.

----

What You Will Learn
-------------------

- How to verify 2D medical image classifiers (single CT/MRI slices)
- How to verify 3D volumetric classifiers using VolumeStar
- How to model clinically relevant perturbations (bias field, gamma, noise)
- The three-stage verification pattern for efficient analysis

Background
----------

Medical imaging models must be robust to real-world imaging variations that
occur during clinical acquisition -- scanner noise, MRI bias field
inhomogeneity, and contrast/brightness changes. Unlike adversarial attacks
in computer vision, these perturbations are **physically motivated** and have
specific mathematical characterizations.

Three-Stage Verification Pattern
----------------------------------

The SPIE 2025 tutorial uses a three-stage approach for efficient verification:

1. **Check original**: Ensure the network correctly classifies the unperturbed image
2. **Bounds check**: Evaluate the network on the lower and upper bound images
3. **Full reachability**: Only run the expensive analysis if stages 1-2 pass

This avoids wasting computation on images that are already misclassified or
trivially non-robust.

2D Classification
------------------

Verify robustness of organ classification from single CT/MRI slices:

.. code-block:: matlab

   % Load network and test image
   load('organ_classifier.mat');
   net = matlab2nnv(trained_net);
   img = test_images(:,:,:,1);
   target = test_labels(1);

   % Stage 1: Check original classification
   y = net.evaluate(single(img));
   [~, pred] = max(y);
   if pred ~= target
       fprintf('Original image misclassified. Skipping verification.\n');
       return;
   end

   % Stage 2: Create perturbation set and check bounds
   epsilon = 0.02;   % L-inf perturbation
   lb = max(single(img) - epsilon, 0);
   ub = min(single(img) + epsilon, 1);

   y_lb = net.evaluate(lb);
   y_ub = net.evaluate(ub);
   [~, pred_lb] = max(y_lb);
   [~, pred_ub] = max(y_ub);

   if pred_lb ~= target || pred_ub ~= target
       fprintf('NOT ROBUST: Classification changes at perturbation bounds.\n');
       return;
   end

   % Stage 3: Full reachability analysis
   IS = ImageStar(lb, ub);
   reachOptions.reachMethod = 'approx-star';
   result = net.verify_robustness(IS, reachOptions, target);

   switch result
       case 1, fprintf('VERIFIED ROBUST under epsilon=%.3f.\n', epsilon);
       case 0, fprintf('NOT ROBUST.\n');
       case 2, fprintf('UNKNOWN (try smaller epsilon or exact-star).\n');
   end

3D Volumetric Classification
------------------------------

For 3D medical images (MRI volumes, CT scans), use VolumeStar:

.. code-block:: matlab

   % Load 3D network and volume
   load('volume_classifier.mat');
   net = matlab2nnv(trained_net);
   volume = test_volumes(:,:,:,:,1);   % H x W x D x C

   % Create VolumeStar (4D perturbation set)
   epsilon = 0.01;
   lb_vol = max(single(volume) - epsilon, 0);
   ub_vol = min(single(volume) + epsilon, 1);
   VS = VolumeStar(lb_vol, ub_vol);

   % Verify
   reachOptions.reachMethod = 'approx-star';
   output_sets = net.reach(VS, reachOptions);
   [out_lb, out_ub] = output_sets.getRanges();

VolumeStar extends ImageStar to 4 dimensions (H x W x D x C), enabling
verification of Conv3D and AveragePooling3D layers used in volumetric models.

Clinically Relevant Perturbations
-----------------------------------

Beyond L-inf noise, medical images are subject to specific acquisition artifacts:

**Bias Field** (MRI):
Intensity inhomogeneity caused by RF coil sensitivity variations. The
perturbation is spatially smooth and multiplicative.

.. code-block:: matlab

   % Bias field perturbation: multiplicative smooth field
   % lb_bias and ub_bias define the range of bias field values per pixel
   IS = ImageStar(img .* lb_bias, img .* ub_bias);

**Gamma Correction** (all modalities):
Brightness/contrast variation from different display settings or scanner
calibration.

.. code-block:: matlab

   % Gamma correction: img^gamma for gamma in [gamma_min, gamma_max]
   lb_gamma = img .^ gamma_max;   % Lower bound (darker)
   ub_gamma = img .^ gamma_min;   % Upper bound (brighter)
   IS = ImageStar(lb_gamma, ub_gamma);

**Random Noise** (all modalities):
Gaussian noise from sensor electronics, especially at high ISO or
low-field MRI.

.. code-block:: matlab

   % Additive noise perturbation
   IS = ImageStar(img - noise_bound, img + noise_bound);

Source Files
^^^^^^^^^^^^

- `Tutorial/SPIE/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/SPIE>`_
- `NN/medmnist/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/medmnist>`_
