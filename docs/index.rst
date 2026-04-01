NNV: Neural Network Verification Toolbox
========================================

.. rst-class:: lead

   **NNV** is a MATLAB toolbox for neural network verification using sound
   set-based reachability analysis. It supports feedforward, convolutional,
   recurrent, graph, and semantic segmentation networks, as well as neural
   ODEs and neural network control systems.

.. rst-class:: lead

   Developed by the `VeriVITAL <https://www.verivital.com/>`_ research group
   at Vanderbilt University.

----

.. grid:: 2 2 3 3
   :gutter: 3

   .. grid-item-card:: Star Set Reachability
      :class-card: sd-border-0 sd-shadow-sm
      :text-align: center

      Exact and approximate reachability analysis with Star sets,
      zonotopes, and abstract domains -- formal safety guarantees.

   .. grid-item-card:: 10+ Architecture Types
      :class-card: sd-border-0 sd-shadow-sm
      :text-align: center

      FFNNs, CNNs, RNNs, SSNNs, Neural ODEs, GNNs, Vision Transformers,
      BNNs, and Neural Network Control Systems.

   .. grid-item-card:: 48 Layer Types
      :class-card: sd-border-0 sd-shadow-sm
      :text-align: center

      Conv1D/2D/3D, ReLU, Sigmoid, Tanh, LSTM, GCN, GINE, BatchNorm,
      MaxPool, and many more -- with exact and approximate methods.

   .. grid-item-card:: ONNX & VNNLIB
      :class-card: sd-border-0 sd-shadow-sm
      :text-align: center

      Load ONNX models, parse VNNLIB specifications, and run
      VNN-COMP benchmarks out of the box.

   .. grid-item-card:: Control Systems
      :class-card: sd-border-0 sd-shadow-sm
      :text-align: center

      Verify closed-loop systems with NN controllers -- linear, nonlinear,
      discrete, and hybrid automaton plant models via CORA integration.

   .. grid-item-card:: VNN-COMP & ARCH-COMP
      :class-card: sd-border-0 sd-shadow-sm
      :text-align: center

      Proven across 11 competition entries in VNN-COMP and ARCH-COMP,
      plus 30+ peer-reviewed publications.

----

Quick Example
-------------

.. code-block:: matlab

   % Load a neural network from ONNX
   dlnet = importNetworkFromONNX('model.onnx', InputDataFormats='BC');
   net = matlab2nnv(dlnet);

   % Define input set as a Star (L-inf ball around a point)
   lb = [0.4; 0.4; 0.4];   % lower bounds
   ub = [0.6; 0.6; 0.6];   % upper bounds
   input_set = Star(lb, ub);

   % Compute reachable output set
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';
   output_sets = net.reach(input_set, reachOptions);

   % Check safety property: output(1) >= 0
   unsafe_region = HalfSpace([1 0], 0);  % y1 <= 0
   result = verify_specification(output_sets, unsafe_region);
   % result: 1 = safe, 0 = unsafe, 2 = unknown

----

.. toctree::
   :maxdepth: 2
   :hidden:

   getting-started/index
   user-guide/index
   nnv3
   application-domains
   examples/index
   theory/index
   publications/index
   testing/index

Getting Started
---------------

.. rst-class:: lead

   Install NNV and run your first neural network verification in minutes.

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: Installation
      :link: getting-started/installation
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Set up NNV with MATLAB toolboxes, CORA, and optional Python
      dependencies for conformal prediction.

   .. grid-item-card:: Quick Start
      :link: getting-started/quickstart
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Create your first input set, run reachability analysis, and check
      safety properties.

   .. grid-item-card:: Courses
      :link: getting-started/courses
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      CS 8395 -- Machine Learning Verification course materials
      from Vanderbilt University.

----

User Guide
----------

.. rst-class:: lead

   Comprehensive guides for all NNV features, from set representations
   to control system verification.

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: Architectures & Layers
      :link: user-guide/architectures
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      All supported network types, 48 layer types, the NN and GNN classes,
      and model loading utilities.

   .. grid-item-card:: Set Representations
      :link: user-guide/set-representations
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Star, ImageStar, VolumeStar, GraphStar, Zono, Box, and more --
      mathematical definitions, constructors, and usage guidance.

   .. grid-item-card:: Verification Methods
      :link: user-guide/verification-methods
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Exact, approximate, probabilistic, and abstract domain methods --
      reachOptions, falsification, and result interpretation.

   .. grid-item-card:: Control Systems
      :link: user-guide/nncs
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Linear, nonlinear, discrete, and hybrid plant models with
      neural network controllers and CORA integration.

   .. grid-item-card:: ONNX & VNNLIB
      :link: user-guide/onnx-vnnlib
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Load ONNX models, parse VNNLIB specs, and convert between
      MATLAB and NNV network formats.

   .. grid-item-card:: LP Solvers
      :link: user-guide/lp-solvers
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      linprog, GLPK, and Gurobi -- solver comparison, selection,
      and performance configuration.

   .. grid-item-card:: Conformal Prediction
      :link: user-guide/conformal-prediction
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Python setup, surrogate model training, and probabilistic
      verification with coverage guarantees.

   .. grid-item-card:: Application Domains
      :link: application-domains
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Aerospace, automotive, medical imaging, power systems,
      cybersecurity, fairness, and more.

   .. grid-item-card:: Docker
      :link: user-guide/docker
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Run NNV in a Docker container with MATLAB.

----

NNV 3.0
-------

.. rst-class:: lead

   See :doc:`what's new in NNV 3.0 <nnv3>` -- VolumeStar, FairNNV,
   probabilistic verification, weight perturbation analysis, and more.

----

Examples & Tutorials
--------------------

.. rst-class:: lead

   Worked examples demonstrating NNV's verification capabilities across
   10+ application domains.

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: Image Classification
      :link: examples/image-classification
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      MNIST, GTSRB, MedMNIST, CIFAR-10 -- adversarial robustness
      verification for image classifiers.

   .. grid-item-card:: ACAS Xu
      :link: examples/acasxu
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Aviation collision avoidance -- verifying all 45 ACAS Xu networks
      with ONNX and VNNLIB.

   .. grid-item-card:: Control Systems
      :link: examples/nncs
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      ACC, AEBS, Inverted Pendulum, DC-DC Buck -- closed-loop safety
      verification with NN controllers.

   .. grid-item-card:: Graph Neural Networks
      :link: examples/gnn
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      GCN and GINE on IEEE bus systems -- verifying neural power
      flow predictions.

   .. grid-item-card:: More Examples
      :link: examples/index
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Semantic segmentation, Neural ODEs, RNNs, malware detection,
      medical imaging, and conference tutorials.

----

Theoretical Foundations
-----------------------

.. rst-class:: lead

   Mathematical details behind NNV's verification algorithms,
   set representations, and probabilistic guarantees.

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Star Set Reachability
      :link: theory/star-reachability
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Star set definitions, exact and approximate reachability algorithms,
      and ReLU splitting strategies.

   .. grid-item-card:: ImageStar & VolumeStar
      :link: theory/imagestar-volumestar
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Extending Star sets to 2D images, 3D volumes, and video data.

   .. grid-item-card:: Probabilistic Verification
      :link: theory/probabilistic
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Conformal inference, surrogate models, and the coverage-confidence
      guarantee framework.

   .. grid-item-card:: Fairness Verification
      :link: theory/fairness
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Counterfactual fairness, individual fairness, and the Verified
      Fairness score.

----

References
----------

The methods implemented in NNV are based upon or used in the following publications:

.. admonition:: Key Publications

   D.\  Manzanas Lopez, S.W. Choi, H.-D. Tran, T.T. Johnson,
   "NNV 2.0: The Neural Network Verification Tool,"
   in *Computer Aided Verification (CAV)*, Springer, 2023.
   `DOI: 10.1007/978-3-031-37703-7_19 <https://doi.org/10.1007/978-3-031-37703-7_19>`__

   H.-D.\  Tran et al.,
   "NNV: A Tool for Verification of Deep Neural Networks and Learning-Enabled
   Autonomous Cyber-Physical Systems,"
   in *Computer Aided Verification (CAV)*, 2020.
   `DOI: 10.1007/978-3-030-53288-8_1 <https://doi.org/10.1007/978-3-030-53288-8_1>`__

   H.-D.\  Tran, S. Bak, W. Xiang, T.T. Johnson,
   "Towards Verification of Large Convolutional Neural Networks Using ImageStars,"
   in *Computer Aided Verification (CAV)*, 2020.

   See the full :doc:`publications list <publications/papers>` and
   :doc:`how to cite <publications/citing>`.

Acknowledgements
----------------

This work is supported in part by AFOSR, DARPA, NSF.
