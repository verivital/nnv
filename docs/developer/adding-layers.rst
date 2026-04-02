Adding New Layer Types
======================

.. rst-class:: lead

   How to implement a custom layer and integrate it into NNV's
   reachability pipeline, based on the actual FullyConnectedLayer
   implementation.

----

The Layer Interface
-------------------

Every NNV layer is a MATLAB ``handle`` class with three key methods:

1. **Constructor**: Store weights, biases, and hyperparameters
2. **evaluate(x)**: Forward pass on a concrete input
3. **reach(varargin)**: Reachability analysis with set-based inputs

The ``reach()`` method receives its arguments via ``varargin`` because
different layers accept different numbers of arguments. The standard
6-argument call from ``NN.reach()`` is:

.. code-block:: matlab

   outputSet = layer.reach(inputSet, method, option, relaxFactor, dis_opt, lp_solver)

However, **graph layers** receive extra arguments (adjacency matrix, edge
features) and **affine layers** ignore the last three (relaxFactor, dis_opt,
lp_solver) since their reachability is exact.

Walkthrough: FullyConnectedLayer
---------------------------------

``FullyConnectedLayer.m`` (~650 lines) is the canonical reference for
implementing a new layer. Here are its key components:

**Properties:**

.. code-block:: matlab

   classdef FullyConnectedLayer < handle
       properties
           Name = 'fully_connected_layer';
           InputSize = 0;
           OutputSize = 0;
           Weights = [];          % W matrix (OutputSize x InputSize)
           Bias = [];             % b vector (OutputSize x 1)
           weightPerturb = [];    % ModelStar weight perturbation spec
       end

**Constructor** (supports 0, 2, or 3 arguments):

.. code-block:: matlab

   function obj = FullyConnectedLayer(varargin)
       switch nargin
           case 2    % FullyConnectedLayer(W, b)
               W = varargin{1}; b = varargin{2};
               obj.Weights = W; obj.Bias = b;
               obj.InputSize = size(W, 2);
               obj.OutputSize = size(W, 1);
           case 3    % FullyConnectedLayer(name, W, b)
               obj.Name = varargin{1};
               W = varargin{2}; b = varargin{3};
               obj.Weights = W; obj.Bias = b;
               obj.InputSize = size(W, 2);
               obj.OutputSize = size(W, 1);
       end
   end

**evaluate(x)** -- flattens multi-dimensional input, then applies y = Wx + b:

.. code-block:: matlab

   function y = evaluate(obj, x)
       n = size(x);
       if length(n) == 2
           x = reshape(x, [n(1)*n(2) 1]);
       elseif length(n) == 3
           x = reshape(x, [n(1)*n(2)*n(3) 1]);
       end
       y = obj.Weights * x + obj.Bias;
   end

**reach(varargin)** -- parses arguments and dispatches:

.. code-block:: matlab

   function IS = reach(varargin)
       % Parse varargin (handles 3 to 7 arguments)
       obj = varargin{1};
       in_images = varargin{2};
       method = varargin{3};
       option = varargin{4};  % 'single' or 'parallel'
       % relaxFactor, dis_opt, lp_solver are accepted but NOT USED
       % (affine reachability is exact, no LP needed)

       if strcmp(method, 'approx-star') || strcmp(method, 'exact-star') ...
               || strcmp(method, 'abs-dom') || contains(method, 'relax-star')
           IS = obj.reach_star_multipleInputs(in_images, option);
       elseif strcmp(method, 'approx-zono')
           IS = obj.reach_zono_multipleInputs(in_images, option);
       end
   end

The Two-Level Reach Pattern
----------------------------

NNV layers follow a standard pattern with two internal methods:

1. ``reach_star_single_input(obj, inputSet)`` -- handles one Star/ImageStar
2. ``reach_star_multipleInputs(obj, inputs, option)`` -- loops over array, handles parallelism

.. code-block:: matlab

   function S = reach_star_multipleInputs(obj, inputs, option)
       if ~iscell(inputs)
           S = obj.reach_star_single_input(inputs);
       else
           n = length(inputs);
           S = cell(1, n);
           if strcmp(option, 'parallel')
               parfor i = 1:n
                   S{i} = obj.reach_star_single_input(inputs{i});
               end
           else
               for i = 1:n
                   S{i} = obj.reach_star_single_input(inputs{i});
               end
           end
       end
   end

For **affine layers**, ``reach_star_single_input`` applies the linear map to
the center and all generators of the Star/ImageStar:

.. code-block:: matlab

   function image = reach_star_single_input(obj, in_image)
       if isa(in_image, 'ImageStar')
           % Apply W to every generator column: V_out = W * V_in
           N = in_image.height * in_image.width * in_image.numChannel;
           n = in_image.numPred;
           V(1,1,:,:) = obj.Weights * reshape(in_image.V, N, n+1);
           V(1,1,:,1) = reshape(V(1,1,:,1), obj.OutputSize, 1) + obj.Bias;
           % Constraints unchanged (affine map preserves constraints)
           image = ImageStar(V, in_image.C, in_image.d, ...
                             in_image.pred_lb, in_image.pred_ub);
       elseif isa(in_image, 'Star')
           image = in_image.affineMap(obj.Weights, obj.Bias);
       end
   end

For **nonlinear layers** (ReLU), ``reach_star_single_input`` is more complex:
it handles VolumeStar→Star→ReLU→VolumeStar conversion, ImageStar→Star→ReLU→ImageStar
conversion, and direct Star input. The actual ReLU logic is in ``PosLin.reach()``.

Implementing a New Layer
--------------------------

To add a new layer type (e.g., a custom activation or normalization):

**Step 1: Create the class file** in ``code/nnv/engine/nn/layers/``:

.. code-block:: matlab

   classdef MyNewLayer < handle
       properties
           Name = 'my_new_layer';
           InputSize = [];
           OutputSize = [];
           % Any learnable parameters
       end

       methods
           function obj = MyNewLayer(varargin)
               % Parse constructor arguments
           end

           function y = evaluate(obj, x)
               % Forward pass on concrete input
           end

           function IS = reach(varargin)
               obj = varargin{1};
               in_sets = varargin{2};
               method = varargin{3};
               option = [];
               if nargin >= 4, option = varargin{4}; end
               % For nonlinear layers, also extract:
               % relaxFactor = varargin{5};
               % dis_opt = varargin{6};
               % lp_solver = varargin{7};

               IS = obj.reach_star_multipleInputs(in_sets, option);
           end
       end
   end

**Step 2: Implement reach_star_single_input**:

- For **affine** operations: apply your transformation to the Star's basis
  matrix ``V`` (center + generators). Constraints ``C``, ``d`` stay unchanged.
- For **nonlinear** operations: you must over-approximate. See how
  ``PosLin.m`` handles ReLU crossing neurons with LP-based relaxation.

**Step 3: Handle multiple set types**:

Your layer may receive Star, ImageStar, VolumeStar, or GraphStar. Common
pattern: convert to Star (``input.toStar()``), apply your logic, convert back.

**Step 4: Register for ONNX import** (optional):

In ``code/nnv/engine/utils/matlab2nnv.m``, add an ``elseif`` clause in the
layer recognition chain:

.. code-block:: matlab

   elseif isa(L, 'nnet.cnn.layer.YourMATLABLayerClass')
       Li = MyNewLayer.parse(L);

The ``parse()`` static method should extract weights/parameters from the
MATLAB layer object and return a new NNV layer instance.

Weight Perturbation Support
-----------------------------

The ``FullyConnectedLayer`` also demonstrates **ModelStar** weight perturbation
support via the ``weightPerturb`` property. When set, the ``reach_star_single_input``
method augments the Star set with additional generators representing the weight
uncertainty. This is an optional extension -- your layer does not need to support
weight perturbation unless you want to analyze parameter uncertainty.

Reference Files
----------------

- ``FullyConnectedLayer.m`` -- canonical affine layer (study this first)
- ``ReluLayer.m`` -- canonical nonlinear layer with exact/approx dispatch
- ``Conv2DLayer.m`` -- convolution with ImageStar preservation
- ``GCNLayer.m`` -- graph-aware layer with extra adjacency argument
- ``PosLin.m`` (in ``funcs/``) -- the actual ReLU math for Star sets
