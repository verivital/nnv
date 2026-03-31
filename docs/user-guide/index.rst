User Guide
==========

.. rst-class:: lead

   Comprehensive guides for all NNV features, from set representations
   to control system verification.

----

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: Architectures & Layers
      :link: architectures
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      All supported network types, 48 layer types, the NN and GNN classes,
      and model loading utilities.

   .. grid-item-card:: Set Representations
      :link: set-representations
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Star, ImageStar, VolumeStar, GraphStar, Zono, Box, and more --
      mathematical definitions, constructors, and usage.

   .. grid-item-card:: Verification Methods
      :link: verification-methods
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Exact, approximate, probabilistic, and abstract domain methods --
      choose the right approach for your problem.

   .. grid-item-card:: Control Systems
      :link: nncs
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Linear, nonlinear, discrete, and hybrid plant models with
      neural network controllers.

   .. grid-item-card:: ONNX & VNNLIB
      :link: onnx-vnnlib
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Load ONNX models, parse VNNLIB specifications, and convert
      between MATLAB and NNV formats.

   .. grid-item-card:: LP Solvers
      :link: lp-solvers
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      linprog, GLPK, and Gurobi -- solver comparison and performance.

   .. grid-item-card:: Conformal Prediction
      :link: conformal-prediction
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Python setup, surrogate model training, and probabilistic
      verification with coverage guarantees.

   .. grid-item-card:: Docker
      :link: docker
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Run NNV in a Docker container with MATLAB.

.. toctree::
   :maxdepth: 2
   :hidden:

   architectures
   set-representations
   verification-methods
   nncs
   onnx-vnnlib
   lp-solvers
   conformal-prediction
   docker
