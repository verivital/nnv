API Reference
=============

.. rst-class:: lead

   Function and class reference for the NNV MATLAB toolbox. Each page
   documents constructor signatures, properties, methods, and usage examples.

----

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: NN Class
      :link: neural-network
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      The main neural network class: reach, evaluate, verify_robustness,
      verify_vnnlib, classify, falsify.

   .. grid-item-card:: GNN Class
      :link: gnn
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Graph neural network class: reach, evaluate, setGraph for
      GCN and GINE architectures.

   .. grid-item-card:: Set Types
      :link: sets
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Star, ImageStar, VolumeStar, GraphStar, Zono, Box, HalfSpace --
      constructors, methods, and properties.

   .. grid-item-card:: Layer Types
      :link: layers
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      All 48 layer types: constructors, evaluate, and reach signatures
      grouped by category.

   .. grid-item-card:: Control Systems
      :link: nncs
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      LinearODE, NonLinearODE, LinearNNCS, NonlinearNNCS --
      plant models and closed-loop composition.

   .. grid-item-card:: Utilities
      :link: utilities
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      matlab2nnv, onnx2nnv, load_vnnlib, verify_specification,
      lpsolver, check_nnv_setup, and more.

   .. grid-item-card:: reachOptions
      :link: reach-options
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Complete reference for the reachOptions struct: all fields,
      types, defaults, and descriptions.

.. toctree::
   :maxdepth: 2
   :hidden:

   neural-network
   gnn
   sets
   layers
   nncs
   utilities
   reach-options
