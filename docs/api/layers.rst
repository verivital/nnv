Layer Types
===========

.. rst-class:: lead

   All 48 layer types supported by NNV. Every layer implements ``evaluate(x)``
   for forward pass and ``reach(inputSet, method, ...)`` for reachability.

----

Common Interface
----------------

All layers follow this interface:

.. code-block:: matlab

   % Forward pass
   y = layer.evaluate(x)

   % Reachability (called internally by NN.reach)
   outputSet = layer.reach(inputSet, method, option, relaxFactor, dis_opt, lp_solver)

   % Static: parse from MATLAB layer (used by matlab2nnv)
   L = LayerClass.parse(matlab_layer)

Input Layers
------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Layer
     - Constructor
     - Notes
   * - ``FeatureInputLayer``
     - ``FeatureInputLayer(inputSize)``
     - Feature vector input
   * - ``ImageInputLayer``
     - ``ImageInputLayer(inputSize)``
     - 2D image input (H x W x C)
   * - ``Image3DInputLayer``
     - ``Image3DInputLayer(inputSize)``
     - 3D volume input (H x W x D x C)
   * - ``SequenceInputLayer``
     - ``SequenceInputLayer(inputSize)``
     - Sequential/time-series input

Linear Layers
-------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Layer
     - Constructor
     - Notes
   * - ``FullyConnectedLayer``
     - ``FullyConnectedLayer(W, b)``
     - Dense: y = Wx + b
   * - ``ElementwiseAffineLayer``
     - ``ElementwiseAffineLayer(scale, offset)``
     - Per-element: y = scale * x + offset

Convolutional Layers
--------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Layer
     - Constructor
     - Notes
   * - ``Conv1DLayer``
     - ``Conv1DLayer(W, b, ...)``
     - 1D convolution with stride, padding, dilation
   * - ``Conv2DLayer``
     - ``Conv2DLayer(W, b, ...)``
     - 2D convolution (workhorse for CNN verification)
   * - ``Conv3DLayer``
     - ``Conv3DLayer(W, b, ...)``
     - 3D convolution for volumetric data
   * - ``TransposedConv1DLayer``
     - ``TransposedConv1DLayer(W, b, ...)``
     - 1D deconvolution
   * - ``TransposedConv2DLayer``
     - ``TransposedConv2DLayer(W, b, ...)``
     - 2D deconvolution (segmentation decoders)

Pooling Layers
--------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Layer
     - Constructor
     - Notes
   * - ``MaxPooling2DLayer``
     - ``MaxPooling2DLayer(poolSize, stride, padding)``
     - 2D max pooling
   * - ``AveragePooling2DLayer``
     - ``AveragePooling2DLayer(poolSize, stride, padding)``
     - 2D average pooling
   * - ``AveragePooling3DLayer``
     - ``AveragePooling3DLayer(...)``
     - 3D average pooling
   * - ``GlobalAveragePooling1DLayer``
     - ``GlobalAveragePooling1DLayer()``
     - Global average over 1D spatial dim
   * - ``GlobalAveragePooling2DLayer``
     - ``GlobalAveragePooling2DLayer()``
     - Global average over 2D spatial dims
   * - ``GlobalAveragePooling3DLayer``
     - ``GlobalAveragePooling3DLayer()``
     - Global average over 3D spatial dims
   * - ``MaxUnpooling2DLayer``
     - ``MaxUnpooling2DLayer(...)``
     - Unpooling for segmentation

Activation Layers
-----------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Layer
     - Constructor
     - Notes
   * - ``ReluLayer``
     - ``ReluLayer()``
     - ReLU: exact-star and approx-star reachability
   * - ``LeakyReluLayer``
     - ``LeakyReluLayer(alpha)``
     - Leaky ReLU: y = max(alpha*x, x)
   * - ``SigmoidLayer``
     - ``SigmoidLayer()``
     - Sigmoid activation
   * - ``TanhLayer``
     - ``TanhLayer()``
     - Hyperbolic tangent
   * - ``HardSigmoidLayer``
     - ``HardSigmoidLayer()``
     - Piecewise-linear sigmoid approximation
   * - ``SaturatingLinearLayer``
     - ``SaturatingLinearLayer()``
     - Clip to [0, 1]
   * - ``SaturatingLinearSymmLayer``
     - ``SaturatingLinearSymmLayer()``
     - Clip to [-1, 1]
   * - ``SignLayer``
     - ``SignLayer()``
     - Sign/step function
   * - ``SoftmaxLayer``
     - ``SoftmaxLayer()``
     - Softmax for classification output

Normalization Layers
--------------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Layer
     - Constructor
     - Notes
   * - ``BatchNormalizationLayer``
     - ``BatchNormalizationLayer(...)``
     - Fused with preceding linear layer during verification
   * - ``LayerNormalizationLayer``
     - ``LayerNormalizationLayer(...)``
     - Layer normalization

Graph Layers
------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Layer
     - Constructor
     - Notes
   * - ``GCNLayer``
     - ``GCNLayer(W, b)``
     - Graph Convolutional Network layer
   * - ``GINELayer``
     - ``GINELayer(W_e, b_e, W1, b1, W2, b2)``
     - Graph Isomorphism Network with Edge features

Recurrent Layers
----------------

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Layer
     - Constructor
     - Notes
   * - ``LstmLayer``
     - ``LstmLayer(Wi, Wh, b, ...)``
     - LSTM cell
   * - ``RecurrentLayer``
     - ``RecurrentLayer(rnn_params)``
     - Simple RNN with Wi, Wh, bh, Wo kernels

Reshaping & Combination Layers
-------------------------------

``FlattenLayer``, ``ReshapeLayer``, ``Resize2DLayer``, ``UpsampleLayer``,
``AdditionLayer``, ``ConcatenationLayer``, ``DepthConcatenationLayer``

Special Layers
--------------

``ODEblockLayer`` (Neural ODE integration), ``PixelClassificationLayer``
(segmentation output), ``LayerS`` (custom layer support)
