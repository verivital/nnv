Developer Guide
===============

.. rst-class:: lead

   Guides for extending NNV to new architectures, set representations,
   and verification problems. Aimed at researchers who want to implement
   their own ideas within the NNV framework.

----

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: Architecture Overview
      :link: architecture
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      How NNV is structured: Computation Engine, Analyzer, directory
      layout, and how NN/GNN/NNCS compose layers.

   .. grid-item-card:: Custom Verification
      :link: custom-verification
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      End-to-end workflows for verifying your own network: FFNN,
      CNN, GNN, NNCS, and custom safety specifications.

   .. grid-item-card:: Adding New Layers
      :link: adding-layers
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      How to implement a new layer type with evaluate and reach
      methods, and register it for ONNX import.

   .. grid-item-card:: Adding New Set Types
      :link: adding-sets
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      How to implement a new set representation with the required
      interface for NNV's reachability pipeline.

   .. grid-item-card:: Adding Reach Methods
      :link: adding-methods
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      How to add a new reachability method and integrate it with
      the existing layer dispatch.

   .. grid-item-card:: Testing & CI/CD
      :link: testing
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Running the test suite, baseline management, code coverage,
      and GitHub Actions CI/CD pipeline.

.. toctree::
   :maxdepth: 2
   :hidden:

   architecture
   custom-verification
   adding-layers
   adding-sets
   adding-methods
   testing
