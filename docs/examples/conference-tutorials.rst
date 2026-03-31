Conference Tutorials
====================

.. rst-class:: lead

   Tutorial materials from major conferences. These provide structured
   introductions to NNV's verification capabilities, suitable for
   workshops and courses.

----

AAAI 2026 -- VNN-COMP Lab
--------------------------

**"The Verification of Neural Networks Competition (VNN-COMP): A Lab for
Benchmark Proposers, Verification Tool Participants, and the Broader AI
Community"** (Lab LH03)

Presented on March 4, 2026 at AAAI-26 in Philadelphia, PA, this half-day
lab introduces participants to neural network verification with hands-on
interactive demos. NNV contributors Taylor T. Johnson and Ben Wooding are
among the organizers. You will learn:

- Introduction to neural network verification and its motivation
- VNN-COMP history, rules, and VNN-COMP 2025 results
- Benchmark formats: VNN-LIB 1.0/2.0 and ONNX
- How to participate as a tool developer or benchmark proposer
- Interactive demos via Google Colab notebooks
- AWS-based evaluation infrastructure overview

- `Slides <https://bit.ly/aaai26vnncomp_slides>`_
- `VNN-LIB Standard <https://www.vnnlib.org/>`_
- `VNN-COMP 2026 <https://github.com/VNN-COMP/vnncomp2026>`_ (co-located with SAIV @ CAV/FLoC 2026, Lisbon, July 24--25)

SPIE 2025 -- Medical Imaging Verification
------------------------------------------

**"Robustness Verification of Medical Imaging Neural Networks"**

A hands-on tutorial covering verification of medical image classifiers and
segmentation networks. You will learn:

- 2D classification robustness with OrganCMNIST
- 3D volumetric classification with VolumeStar
- Segmentation robustness under bias field, gamma correction, and noise
- When to use conformal prediction for large segmentation models

- Location: ``examples/Tutorial/SPIE/``
- See :doc:`medical-imaging` for the worked examples

DSN 2024 -- Dependable Systems
-------------------------------

**"Tutorial: Safe, Secure, and Trustworthy AI via Formal Verification of
Neural Networks and Autonomous CPS with NNV"**

This tutorial covers NNV's core verification workflow with emphasis on
safety-critical systems. You will learn:

- Loading and converting neural networks (MATLAB and ONNX formats)
- Defining input perturbation sets (Star, ImageStar)
- Running exact and approximate reachability analysis
- Verifying closed-loop control systems (NNCS)
- Adversarial robustness analysis for image classifiers

- `Paper: DOI 10.1109/DSN-S60304.2024.00027 <https://doi.org/10.1109/DSN-S60304.2024.00027>`_

EMSOFT 2023 -- Embedded Software
---------------------------------

**"Tutorial: Neural Network and Autonomous CPS Formal Verification for
Trustworthy AI and Safe Autonomy"**

An introduction to NNV focusing on embedded and cyber-physical systems.
You will learn:

- Star set reachability fundamentals
- FFNN and CNN verification
- Neural network control system verification (ACC benchmark)
- Interpreting verification results for safety-critical embedded systems

- `Paper: DOI 10.1145/3607890.3608454 <https://doi.org/10.1145/3607890.3608454>`_
- Location: ``examples/Tutorial/``

Running Tutorials Online
--------------------------

All tutorials can be run without local installation:

- **MATLAB Online**: `Try NNV on MATLAB Online <https://matlab.mathworks.com/>`_
  (a MATLAB license may be required for some examples, but many run as guest)
- **CodeOcean**: `NNV CAV 2023 capsule <https://doi.org/10.24433/CO.0803700.v1>`_
  (runs in browser, no installation needed)

All tutorial source code is available at
`examples/Tutorial/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial>`_.
