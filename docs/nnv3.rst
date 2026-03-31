NNV 3.0 Summary
================

.. rst-class:: lead

   NNV 3.0 introduces major new verification capabilities, building on
   NNV 1.0 (CAV 2020) and NNV 2.0 (CAV 2023). This page summarizes
   what's new and links to the relevant documentation.

----

Version History
---------------

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Version
     - Venue
     - Major Additions
   * - **NNV 1.0**
     - CAV 2020
     - Star-based reachability for FFNNs and CNNs, ImageStar representation,
       ReLU verification, NNCS verification
   * - **NNV 2.0**
     - CAV 2023
     - RNN/LSTM support, semantic segmentation networks, Neural ODEs,
       ONNX and VNNLIB format support, VNN-COMP participation, relaxed reachability
   * - **NNV 3.0**
     - FM 2026
     - VolumeStar, FairNNV, probabilistic verification, weight perturbation
       (ModelStar), time-dependent networks, GNNs, malware detection benchmark

What's New in NNV 3.0
---------------------

VolumeStar: 3D/Video Data Verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First-of-its-kind verification for video classification and 3D medical imaging
networks. The VolumeStar extends the Star/ImageStar hierarchy to 4D tensors
(H x W x C x F), enabling robustness verification of C3D, I3D, and 3D-CNN architectures.

- Handles both spatial AND temporal perturbations
- Supports Conv3D and AveragePooling3D layers
- Applications: video classification, 3D MRI/CT analysis

**Documentation:** :doc:`user-guide/set-representations` (VolumeStar section),
:doc:`theory/imagestar-volumestar`, :doc:`examples/medical-imaging`

FairNNV: Fairness Verification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Formal verification of algorithmic fairness via reachability analysis, proving
that model predictions are invariant (or bounded) with respect to sensitive
attributes like gender or race.

- **Counterfactual fairness**: Predictions unchanged when sensitive attributes are altered
- **Individual fairness**: Similar individuals receive similar predictions
- **Verified Fairness (VF) score**: Quantitative, formally grounded fairness measure
- Guarantees hold for ALL inputs, not just statistical samples

**Documentation:** :doc:`theory/fairness`

Probabilistic Verification via Conformal Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Data-driven reachability based on conformal inference for problems where full
worst-case verification is computationally intractable (e.g., large semantic
segmentation networks).

- User specifies perturbation space, sampling distribution, miscoverage level
- Returns probabilistic output sets with formal coverage guarantees
- Python surrogate model training (linear or ReLU)

**Documentation:** :doc:`user-guide/conformal-prediction`, :doc:`theory/probabilistic`

ModelStar: Weight Perturbation Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Verification under parameter uncertainties caused by quantization, hardware
errors, or device variations. Star set augmentations capture interval-based
parameter perturbations for linear layers.

- Single-layer perturbation: exact Star set representation
- Multi-layer perturbation: sound over-approximation
- Implemented for fully-connected and Conv2D layers

**Documentation:** :doc:`user-guide/verification-methods`

Time-Dependent Neural Networks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Extends RNN verification beyond fixed-length sequences to handle variable-length
time series with time horizon T in [T_min, T_max].

- Computes union of reachable sets over possible time horizons
- Applications: NLP, sensor data analysis, sequential decision systems

**Documentation:** :doc:`user-guide/architectures`, :doc:`examples/rnn`

Graph Neural Networks (GCN, GINE)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New ``GNN`` class and ``GraphStar`` set representation for verifying graph
neural networks under node feature perturbations.

- GCN (Graph Convolutional Network) layers
- GINE (Graph Isomorphism Network with Edge features) layers
- Demonstrated on IEEE bus power systems (24, 39, 118 buses)

**Documentation:** :doc:`user-guide/architectures` (GNN section),
:doc:`user-guide/set-representations` (GraphStar section), :doc:`examples/gnn`

Malware Detection Benchmark
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New cybersecurity verification domain: robustness verification of DNN-based
malware classifiers against adversarial evasion attacks.

- Feature-space and image-based perturbations
- Formal proof that modifications cannot cause misclassification

**Documentation:** :doc:`examples/malware`

Software Engineering Improvements
----------------------------------

- **check_nnv_setup()**: Comprehensive diagnostic tool with three-tier dependency checking
- **NNVVERSION()**: Programmatic version querying
- **CI/CD pipeline**: GitHub Actions with automated testing (833 tests, 100% pass rate)
- **QuickStart folder**: ``test_installation.m`` and ``simple_verification.m``
- **48.6% code coverage** (86.2% in nn/, 88.9% in set/)
- **Regression testing** with baseline comparison (46 baselines)

Comparison vs Other Tools
-------------------------

See the full comparison table in :doc:`user-guide/architectures` for a
14-tool comparison across all supported architectures and applications.

NNV3 is the only tool that supports **all** of the following: GNNs, TDNNs,
3D CNNs/video, weight perturbation analysis, fairness verification, and
probabilistic verification -- in addition to standard FFNN, CNN, RNN, SSNN,
Neural ODE, and NNCS verification

Citation
--------

See :doc:`publications/citing` for BibTeX entries for NNV 1.0, 2.0, and 3.0.
