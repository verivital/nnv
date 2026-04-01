Application Domains
===================

.. rst-class:: lead

   NNV has been applied to verify neural networks across 10+ safety-critical
   domains. This page summarizes each domain, the verification challenges
   involved, and pointers to examples and publications.

----

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: Aerospace
      :class-card: sd-border-0 sd-shadow-sm

      Collision avoidance advisories (ACAS Xu)

   .. grid-item-card:: Automotive
      :class-card: sd-border-0 sd-shadow-sm

      Cruise control, emergency braking, perception

   .. grid-item-card:: Medical Imaging
      :class-card: sd-border-0 sd-shadow-sm

      Organ classification, 3D segmentation, CT/MRI

   .. grid-item-card:: Power Systems
      :class-card: sd-border-0 sd-shadow-sm

      Optimal power flow, voltage safety (GNNV)

   .. grid-item-card:: Cybersecurity
      :class-card: sd-border-0 sd-shadow-sm

      Malware detection, adversarial evasion

   .. grid-item-card:: Fairness / AI Ethics
      :class-card: sd-border-0 sd-shadow-sm

      Algorithmic fairness in credit scoring, hiring

   .. grid-item-card:: Financial Predictions
      :class-card: sd-border-0 sd-shadow-sm

      Credit approval, loan risk, fairness certification

   .. grid-item-card:: Robotics & Control
      :class-card: sd-border-0 sd-shadow-sm

      Inverted pendulum, DC-DC buck converter

   .. grid-item-card:: Industrial Prognostics
      :class-card: sd-border-0 sd-shadow-sm

      Remaining useful life prediction (Collins RUL)

   .. grid-item-card:: Autonomous Driving
      :class-card: sd-border-0 sd-shadow-sm

      Traffic sign recognition, street scene segmentation

----

Aerospace
---------

**Problem:** Unmanned aircraft collision avoidance systems use neural networks
to issue turn advisories based on intruder position and heading. A wrong
advisory could lead to a mid-air collision.

**Networks:** 45 ACAS Xu feed-forward networks (5 inputs, 5 outputs, 6 hidden
layers of 50 neurons each).

**Verification approach:** Load ONNX models and VNNLIB property specifications,
verify safety properties (e.g., "if the intruder is directly ahead, do not
advise Clear-of-Conflict") using ``approx-star`` and ``exact-star``.

**NNV features used:** Star sets, ONNX/VNNLIB support, HalfSpace safety specs

.. list-table::
   :widths: 30 70

   * - Examples
     - :doc:`/examples/acasxu`
   * - Key publication
     - Lopez et al., "Evaluation of Neural Network Verification Methods for Air
       to Air Collision Avoidance", AIAA JAT 2022

----

Automotive
----------

**Problem:** Neural network controllers in vehicles must maintain safe distances
(adaptive cruise control) and brake in emergencies (AEBS). A control failure
could cause a rear-end collision.

**Networks:** Feed-forward controllers with ReLU and saturation activations,
embedded in closed-loop systems with continuous/discrete plant dynamics.

**Verification approach:** Compose NN controller with plant model
(``NonlinearNNCS``, ``DLinearNNCS``), compute reachable state sets over time,
verify safety constraints (e.g., following distance > 10m + 1.4s time gap).

**NNV features used:** NNCS classes, LinearODE/NonLinearODE plants, CORA integration,
temporal reachability

.. list-table::
   :widths: 30 70

   * - Examples
     - :doc:`/examples/nncs` (ACC, AEBS)
   * - Competition
     - ARCH-COMP AINNCS category (2019--2024)

**Traffic sign and scene perception:** GTSRB traffic sign robustness, CamVid
street scene segmentation -- verifying that perception networks are robust
to adversarial perturbations.

.. list-table::
   :widths: 30 70

   * - Examples
     - :doc:`/examples/image-classification` (GTSRB),
       :doc:`/examples/semantic-segmentation` (CamVid)

----

Medical Imaging
---------------

**Problem:** Neural networks for organ classification and tissue segmentation
in CT/MRI scans must be robust to realistic imaging artifacts (noise, bias
fields, contrast variations). A misclassification could lead to incorrect
diagnosis.

**Networks:** 2D CNNs (single slice) and 3D CNNs (volumetric) for
classification and segmentation of medical images.

**Verification approach:** Use ImageStar (2D) and VolumeStar (3D) to model
clinically relevant perturbations -- L-inf noise, MRI bias field
inhomogeneity, gamma correction. For large segmentation networks, conformal
prediction provides probabilistic coverage guarantees.

**NNV features used:** ImageStar, VolumeStar, Conv3D layers, conformal prediction

.. list-table::
   :widths: 30 70

   * - Examples
     - :doc:`/examples/medical-imaging` (SPIE 2025), :doc:`/examples/image-classification` (MedMNIST)
   * - Key publications
     - Hashemi et al., "Scaling Data-Driven Probabilistic Robustness Analysis
       for Semantic Segmentation Neural Networks", NeurIPS 2025;
       Sasaki et al., "Robustness Verification of Video Classification Neural
       Networks", FormaliSE 2025

----

Power Systems
-------------

**Problem:** Graph neural networks predict voltage levels and optimal power
flow in electrical grids. Incorrect predictions under input uncertainty
(load fluctuations, renewable variability) could violate voltage safety
limits and damage equipment.

**Networks:** GCN (Graph Convolutional Network) and GINE (Graph Isomorphism
Network with Edge features) trained on IEEE bus systems.

**Verification approach:** Use GraphStar to model node feature perturbations
(power injection uncertainties), propagate through GNN layers, verify that
all predicted voltages remain within [0.95, 1.05] per unit.

**NNV features used:** GNN class, GraphStar, GCNLayer, GINELayer

.. list-table::
   :widths: 30 70

   * - Examples
     - :doc:`/examples/gnn` (IEEE 24-bus, 39-bus, 118-bus)
   * - Benchmarks
     - Optimal Power Flow, AC Power Flow prediction

----

Cybersecurity
-------------

**Problem:** DNN-based malware classifiers are vulnerable to adversarial
evasion attacks where malware authors make small feature modifications to
evade detection, allowing malicious software to pass as benign.

**Networks:** Feed-forward classifiers trained on the BODMAS malware dataset
with feature-based and image-based input representations.

**Verification approach:** Feature-aware adversarial perturbations -- only
modify features that an attacker can realistically change (e.g., discrete
features with large range). Two-phase verification: random falsification
search, then full Star reachability.

**NNV features used:** Star sets, feature-aware perturbation bounds, ONNX loading

.. list-table::
   :widths: 30 70

   * - Examples
     - :doc:`/examples/malware`
   * - Key publication
     - Robinette et al., "Case Study: Neural Network Malware Detection
       Verification for Feature and Image Datasets", FormaliSE 2024

----

Fairness & AI Ethics
--------------------

**Problem:** AI models used for high-stakes decisions (loan approvals, credit
scoring, hiring) may discriminate based on sensitive attributes like gender
or race. Statistical auditing only tests sampled inputs and can miss
systematic bias.

**Verification approach:** FairNNV uses Star-based reachability to formally
verify fairness properties over continuous input regions. Counterfactual
fairness (predictions unchanged when sensitive attributes are altered) and
individual fairness (similar inputs get similar outputs) are checked
exhaustively, not just on samples.

**NNV features used:** Star sets, FairNNV module, Verified Fairness (VF) score

.. list-table::
   :widths: 30 70

   * - Theory
     - :doc:`/theory/fairness`
   * - Key publication
     - Tumlin et al., "FairNNV: The Neural Network Verification Tool For
       Certifying Fairness", ICAIF 2024

----

Financial Predictions
---------------------

**Problem:** Neural networks used for credit approval, loan risk assessment, and
financial scoring must be both accurate and fair. Regulatory and ethical
requirements demand that predictions do not systematically discriminate based
on protected attributes such as race, gender, or age.

**Networks:** Feed-forward classifiers trained on tabular financial data
(e.g., Adult Census, German Credit, Bank Marketing datasets).

**Verification approach:** FairNNV verifies counterfactual and individual fairness
over continuous input regions, certifying that decisions remain consistent
when sensitive attributes are perturbed. Unlike statistical auditing over finite
test sets, reachability-based analysis provides formal guarantees. Evaluation on
the Adult Census dataset shows 87--89% Verified Fairness scores under
counterfactual specifications, with verification times below 0.03s per sample.

**NNV features used:** Star sets, FairNNV module, Verified Fairness (VF) score

.. list-table::
   :widths: 30 70

   * - Theory
     - :doc:`/theory/fairness`
   * - Key publication
     - Tumlin et al., "FairNNV: The Neural Network Verification Tool For
       Certifying Fairness", ICAIF 2024

----

Robotics & Embedded Control
----------------------------

**Problem:** Neural network controllers in robotic and power electronic systems
must maintain stability and safety under all operating conditions. Instability
in an inverted pendulum or voltage regulator can cause physical damage.

**Networks:** Small feed-forward controllers in feedback loops with linear
or discrete plant models.

**Verification approach:** Compose NN controller with plant dynamics
(``DLinearNNCS``), verify stability and safety over time horizons using
exact-star or approx-star reachability.

**NNV features used:** NNCS classes, DLinearODE, exact-star reachability

.. list-table::
   :widths: 30 70

   * - Examples
     - :doc:`/examples/nncs` (Inverted Pendulum, DC-DC Buck converter)

----

Industrial Prognostics
-----------------------

**Problem:** CNN-based Remaining Useful Life (RUL) predictors estimate when
equipment will fail. A key requirement is **monotonicity** -- the predicted
RUL must decrease over time. Weight perturbations from quantization or
hardware errors could violate this property.

**Verification approach:** Weight perturbation analysis (ModelStar) --
perturb a fraction of network weights and verify that the monotonicity
property is preserved. Find the maximum tolerable perturbation level.

**NNV features used:** Weight perturbation (WPutils), ONNX loading, VNNLIB
monotonicity properties

.. list-table::
   :widths: 30 70

   * - Examples
     - `Tutorial/NN/Collins_RUL/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NN/Collins_RUL>`_
   * - Key publication
     - Zubair et al., "Verification of Neural Network Robustness Against Weight
       Perturbations Using Star Sets", IEEE CAI 2025

----

Autonomous Dynamics Learning
------------------------------

**Problem:** Neural ODEs learn continuous-time dynamics from data. When used
for prediction or control, the learned dynamics must be verified to stay
within safe state regions over a time horizon.

**Networks:** ODEblockLayer with learned linear or nonlinear dynamics,
followed by feed-forward layers.

**Verification approach:** Use ``direct`` (matrix exponential) for linear ODE
blocks or ``ode45`` (numerical integration via CORA) for nonlinear blocks.
Verify that predicted trajectories remain within safe bounds.

**NNV features used:** ODEblockLayer, LinearODE/NonLinearODE, CORA integration

.. list-table::
   :widths: 30 70

   * - Examples
     - :doc:`/examples/neural-odes` (Spiral, CartPole, FPA)
   * - Key publication
     - Lopez et al., "Reachability Analysis of a General Class of Neural
       Ordinary Differential Equations", FORMATS 2022
