Supported Architectures & Layers
================================

.. rst-class:: lead

   NNV supports verification of 10+ neural network architecture types
   through a unified ``NN`` class and 48 individual layer types.

----

Network Types
-------------

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Architecture
     - Description
     - Key Set Type
   * - Feed-Forward (FFNN)
     - Fully-connected layers with activation functions
     - Star
   * - Convolutional (CNN)
     - Conv2D, pooling, batch normalization for image classification
     - ImageStar
   * - Recurrent (RNN/LSTM)
     - Sequential data with recurrent/LSTM cells
     - Star (per time step)
   * - Semantic Segmentation
     - Pixel-level classification (transposed conv, unpooling)
     - ImageStar
   * - Neural ODE
     - Continuous-time dynamics learned as ODE blocks
     - Star + Zono (via CORA)
   * - Graph Neural Network
     - GCN and GINE layers for graph-structured data
     - GraphStar
   * - Binary Neural Network
     - Networks with binary weights/activations
     - Star
   * - Vision Transformer
     - Attention-based image classification
     - Star
   * - 3D / Volumetric CNN
     - Conv3D for medical images and video classification
     - VolumeStar
   * - Time-Dependent NN
     - Variable-length sequence verification
     - Star (union over time horizons)

Feature Overview
^^^^^^^^^^^^^^^^

The table below summarizes the major features available across NNV versions.
Items in **bold** were introduced in NNV3.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Feature
     - Supported (NNV 1.0, NNV 2.0, NNV3)
   * - Network Types
     - FFNN, CNN, NeuralODE, SSNN, RNN, **GNN**, **TDNN**, **3D CNN**
   * - Layers
     - MaxPool, Conv, BN, AvgPool, FC, MaxUnpool, TC, DC, NODE, **GCN**, **GINE**, **Conv3D**
   * - Activations
     - ReLU, Satlin, Sigmoid, Tanh, Leaky ReLU, Satlins
   * - Plant Dynamics
     - Linear ODE, Nonlinear ODE, Continuous & Discrete Time, Hybrid Automata
   * - Set Representations
     - Polyhedron, Zonotope, Star, ImageStar, **VolumeStar**, **ModelStar**, **GraphStar**
   * - Reach Methods
     - exact, approx, abs-dom, relax-*
   * - Visualization
     - Exact and over-approximate reachable sets
   * - Verification
     - Safety, Robustness, VNNLIB, **Fairness**, **Weight Perturbation**, **Probabilistic**
   * - Miscellaneous
     - Parallel computing, counterexample generation, ONNX, **CI/CD**

NNV vs Other Tools
^^^^^^^^^^^^^^^^^^

The following table compares NNV3 against 13 other verification tools across
supported architectures and applications.
A **~** indicates possible or partial support.

.. raw:: html

   <div class="table-scroll">

.. list-table::
   :header-rows: 1
   :class: comparison-table

   * -
     - NNV3
     - a,b-Crown
     - CORA
     - FastBATLLNN
     - Marabou
     - NeuralSAT
     - NeVer2
     - nnenum
     - ReachNN
     - Reluplex
     - SobolBox
     - Verinet
     - Verisig
     - StarV
   * - FFNN
     - Yes
     - Yes
     - Yes
     - Yes
     - Yes
     - Yes
     - Yes
     - Yes
     - \-
     - Yes
     - Yes
     - Yes
     - \-
     - Yes
   * - CNN
     - Yes
     - Yes
     - Yes
     - \-
     - Yes
     - Yes
     - Yes
     - Yes
     - \-
     - \-
     - Yes
     - Yes
     - \-
     - Yes
   * - RNN
     - Yes
     - Yes
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - Yes
   * - SSNN
     - Yes
     - ~
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - Yes
   * - TDNNs
     - Yes
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
   * - GNNs
     - Yes
     - \-
     - Yes
     - \-
     - ~
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
   * - 3D Data
     - Yes
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
   * - Weight Perturbation
     - Yes
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
   * - Neural ODEs
     - Yes
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - \-
     - Yes
   * - NNCS
     - Yes
     - ~
     - Yes
     - \-
     - \-
     - \-
     - \-
     - \-
     - Yes
     - \-
     - \-
     - \-
     - Yes
     - Yes

.. raw:: html

   </div>

The NN Class
------------

The ``NN`` class is the unified interface for all feed-forward-style networks
(FFNNs, CNNs, segmentation nets, BNNs, etc.):

.. code-block:: matlab

   % Create from a cell array of layers
   layers = {FullyConnectedLayer(W1, b1), ReLULayer(), FullyConnectedLayer(W2, b2)};
   net = NN(layers);

   % Or convert from MATLAB/ONNX format
   dlnet = importNetworkFromONNX('model.onnx', InputDataFormats='BC');
   net = matlab2nnv(dlnet);

**Key Properties:**

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Property
     - Type
     - Description
   * - ``Layers``
     - cell array
     - Ordered list of layer objects
   * - ``Connections``
     - table
     - DAG connectivity (for non-sequential architectures)
   * - ``Name``
     - string
     - Network name
   * - ``numLayers``
     - int
     - Number of layers
   * - ``numNeurons``
     - int
     - Total neuron count
   * - ``InputSize``
     - array
     - Input dimensions
   * - ``OutputSize``
     - array
     - Output dimensions

**Key Methods:**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Method
     - Description
   * - ``evaluate(x)``
     - Forward pass on a concrete input vector/image
   * - ``reach(inputSet, reachOptions)``
     - Compute reachable output set(s) from an input set
   * - ``verify(inputSet, property, reachOptions)``
     - Verify a safety property over an input region
   * - ``verify_vnnlib(vnnlib_file, reachOptions)``
     - Verify properties from a VNNLIB specification file

The GNN Class
-------------

The ``GNN`` class handles graph neural networks with adjacency-based message passing:

.. code-block:: matlab

   % Create GNN with layers and graph structure
   gnn = GNN(layers, A_norm, adj_list);

   % Compute reachable set for node feature perturbation
   input_set = GraphStar(NF, LB, UB);
   output_sets = gnn.reach(input_set, reachOptions);

**Key Properties:**

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Property
     - Type
     - Description
   * - ``Layers``
     - cell array
     - GCN/GINE and activation layers
   * - ``A_norm``
     - matrix
     - Normalized adjacency matrix
   * - ``adj_list``
     - matrix
     - Edge list representation
   * - ``E``
     - matrix
     - Edge feature matrix (for GINE layers)
   * - ``edge_weights``
     - vector
     - Optional edge weight vector

Layer Reference
---------------

NNV implements 48 layer types organized into the following categories.

Input Layers
^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``FeatureInputLayer``
     - Input layer for feature vectors
   * - ``ImageInputLayer``
     - Input layer for 2D images (H x W x C)
   * - ``Image3DInputLayer``
     - Input layer for 3D volumes (H x W x D x C)
   * - ``SequenceInputLayer``
     - Input layer for sequential/time-series data
   * - ``PlaceholderLayer``
     - Generic placeholder for unsupported import layers

Linear / Affine Layers
^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``FullyConnectedLayer(W, b)``
     - Dense layer: y = Wx + b
   * - ``ElementwiseAffineLayer``
     - Per-element scale and offset: y = scale * x + offset

Convolutional Layers
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``Conv1DLayer``
     - 1D convolution with stride, padding, dilation
   * - ``Conv2DLayer``
     - 2D convolution for images (the workhorse of CNN verification)
   * - ``Conv3DLayer``
     - 3D convolution for volumetric data / video frames
   * - ``TransposedConv1DLayer``
     - 1D transposed convolution (deconvolution)
   * - ``TransposedConv2DLayer``
     - 2D transposed convolution (used in segmentation decoders)

Pooling Layers
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``MaxPooling2DLayer``
     - 2D max pooling with configurable pool size and stride
   * - ``AveragePooling2DLayer``
     - 2D average pooling
   * - ``AveragePooling3DLayer``
     - 3D average pooling for volumetric data
   * - ``GlobalAveragePooling1DLayer``
     - Global average over 1D spatial dimension
   * - ``GlobalAveragePooling2DLayer``
     - Global average over 2D spatial dimensions
   * - ``GlobalAveragePooling3DLayer``
     - Global average over 3D spatial dimensions
   * - ``MaxUnpooling2DLayer``
     - Unpooling for segmentation decoder networks

Activation Layers
^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``ReluLayer``
     - ReLU: y = max(0, x) -- supports exact and approximate-star reachability
   * - ``LeakyReluLayer``
     - Leaky ReLU: y = max(alpha*x, x)
   * - ``SigmoidLayer``
     - Sigmoid: y = 1/(1 + exp(-x))
   * - ``TanhLayer``
     - Hyperbolic tangent
   * - ``HardSigmoidLayer``
     - Piecewise-linear approximation of sigmoid
   * - ``SaturatingLinearLayer``
     - Clipped linear: y = clip(x, 0, 1)
   * - ``SaturatingLinearSymmLayer``
     - Symmetric saturation: y = clip(x, -1, 1)
   * - ``SignLayer``
     - Sign/step function
   * - ``SoftmaxLayer``
     - Softmax for classification output
   * - ``ActivationFunctionLayer``
     - Base class for custom activation functions

Normalization Layers
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``BatchNormalizationLayer``
     - Batch normalization (fused with preceding conv/FC during verification)
   * - ``LayerNormalizationLayer``
     - Layer normalization

Graph Neural Network Layers
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``GCNLayer``
     - Graph Convolutional Network layer (Kipf & Welling)
   * - ``GINELayer``
     - Graph Isomorphism Network with Edge features

Recurrent Layers
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``LstmLayer``
     - Long Short-Term Memory cell with input/recurrent/output gates
   * - ``RecurrentLayer``
     - Generic simple recurrent layer (Wi, Wh kernels)

Reshaping & Structural Layers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``FlattenLayer``
     - Flatten multi-dimensional input to 1D vector
   * - ``ReshapeLayer``
     - Reshape tensor dimensions
   * - ``ReshapeToConcatenationLayer``
     - Reshape for concatenation compatibility
   * - ``Resize2DLayer``
     - Resize spatial dimensions (interpolation)
   * - ``UpsampleLayer``
     - Upsample spatial dimensions

Combination Layers
^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``AdditionLayer``
     - Element-wise addition of multiple inputs (for skip connections)
   * - ``ConcatenationLayer``
     - Channel/feature concatenation
   * - ``DepthConcatenationLayer``
     - Depth-wise concatenation

Special Layers
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Layer
     - Description
   * - ``ODEblockLayer``
     - Neural ODE integration block (wraps LinearODE/NonLinearODE)
   * - ``PixelClassificationLayer``
     - Pixel-level classification output (for segmentation)
   * - ``LayerS``
     - Custom/generic layer support

Model Loading Utilities
-----------------------

NNV provides several utilities for loading models from different formats:

.. code-block:: matlab

   % From ONNX format
   dlnet = importNetworkFromONNX('model.onnx', InputDataFormats='BC');
   net = matlab2nnv(dlnet);

   % From MATLAB deep learning format
   net = matlab2nnv(trained_network);  % SeriesNetwork, DAGNetwork, or dlnetwork

   % Using the dedicated ONNX converter
   net = onnx2nnv('model.onnx');

   % From .mat file
   net = load_NN_from_mat('network.mat');

See :doc:`onnx-vnnlib` for detailed ONNX and VNNLIB workflow documentation.
