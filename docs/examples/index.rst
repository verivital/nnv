Examples & Tutorials
====================

.. rst-class:: lead

   Worked examples demonstrating NNV's verification capabilities across
   10+ application domains, from image classification to power systems.

----

Recommended Learning Path
--------------------------

1. **Start here**: :doc:`/getting-started/quickstart` -- your first verification
2. **Image classification**: :doc:`image-classification` -- MNIST, GTSRB robustness
3. **Safety-critical systems**: :doc:`acasxu` -- ACAS Xu with ONNX/VNNLIB
4. **Control systems**: :doc:`nncs` -- ACC, AEBS closed-loop verification
5. **Advanced topics**: GNNs, Neural ODEs, probabilistic verification

.. grid:: 1 2 2 3
   :gutter: 3

   .. grid-item-card:: Image Classification
      :link: image-classification
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      MNIST, GTSRB, MedMNIST, CIFAR-10 -- adversarial robustness
      verification for image classifiers.

   .. grid-item-card:: ACAS Xu
      :link: acasxu
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Aviation collision avoidance -- verifying all 45 ACAS Xu networks
      with ONNX and VNNLIB.

   .. grid-item-card:: Semantic Segmentation
      :link: semantic-segmentation
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Pixel-level verification for M2NIST and CamVid segmentation
      networks.

   .. grid-item-card:: Graph Neural Networks
      :link: gnn
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      GCN and GINE on IEEE bus systems for power grid safety
      verification.

   .. grid-item-card:: Neural ODEs
      :link: neural-odes
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Reachability of continuous-time dynamics: spiral, CartPole,
      FPA systems.

   .. grid-item-card:: Recurrent Neural Networks
      :link: rnn
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Sequence-based robustness verification with variable
      time horizons.

   .. grid-item-card:: Control Systems
      :link: nncs
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      ACC, AEBS, Inverted Pendulum, DC-DC Buck converter --
      closed-loop safety verification.

   .. grid-item-card:: Malware Detection
      :link: malware
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      Adversarial robustness of DNN-based malware classifiers
      (BODMAS dataset).

   .. grid-item-card:: Medical Imaging
      :link: medical-imaging
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      SPIE 2025 tutorial: 2D/3D classification and segmentation
      robustness for medical images.

   .. grid-item-card:: Conference Tutorials
      :link: conference-tutorials
      :link-type: doc
      :class-card: sd-border-0 sd-shadow-sm

      SPIE 2025, DSN 2024, EMSOFT 2023 -- tutorial materials
      and references.

.. toctree::
   :maxdepth: 2
   :hidden:

   image-classification
   acasxu
   semantic-segmentation
   gnn
   neural-odes
   rnn
   nncs
   malware
   medical-imaging
   conference-tutorials
