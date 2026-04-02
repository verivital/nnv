Adding New Set Types
====================

.. rst-class:: lead

   How to implement a new set representation and integrate it with
   NNV's reachability pipeline.

----

How Sets Work in NNV
---------------------

NNV's set types all follow the same core idea: represent a region of possible
values using a **center point**, **generator directions** (basis vectors), and
**linear constraints** on how far you can move in each direction.

The fundamental representation is:

.. math::

   S = \{ x \mid x = V(:,1) + \sum_{i=1}^{m} \alpha_i \cdot V(:,i+1), \;\; C\alpha \leq d \}

where ``V(:,1)`` is the center, ``V(:,2:end)`` are generators, and ``C``, ``d``
define the feasible region for the predicate variables ``alpha``.

Each specialized set type (ImageStar, VolumeStar, GraphStar) stores ``V`` in a
different shape to preserve the data structure (images, volumes, graphs) while
using the same constraint representation.

Required Methods
-----------------

To work with NNV's reachability pipeline, a set type needs:

**Core (required by all layers):**

.. code-block:: matlab

   s2 = s.affineMap(W, b)      % Linear transformation: y = Wx + b
   [lb, ub] = s.getRanges()    % Element-wise bounds (via LP)
   bool = s.isEmptySet()       % Emptiness check
   s_star = s.toStar()         % Convert to flat Star (for activation layers)

**For verification:**

.. code-block:: matlab

   bool = s.contains(point)    % Point membership test
   points = s.sample(N)        % Random sampling

**For usability:**

.. code-block:: matlab

   s.plot()                    % Visualization

The ``toStar()`` Method
^^^^^^^^^^^^^^^^^^^^^^^^

This is the most critical method for integration. Activation layers (ReLU,
Sigmoid, etc.) operate on flat vectors via the ``PosLin`` class, which expects
a Star set. Your custom set must be convertible to Star and back:

.. code-block:: matlab

   % In your set class:
   function s = toStar(obj)
       % Reshape V from (domain-specific shape) to (flat_dim x numBasis)
       flat_V = reshape(obj.V, [], size(obj.V, ndims(obj.V)));
       s = Star(flat_V, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
   end

   % Static method to convert back:
   function gs = fromStar(star, dim1, dim2, ...)
       new_V = reshape(star.V, [dim1, dim2, ..., star.nVar + 1]);
       gs = MyCustomSet(new_V, star.C, star.d, ...
                        star.predicate_lb, star.predicate_ub);
   end

This is exactly how ``GraphStar`` works in ``GNN.reach()``:

.. code-block:: matlab

   % Before ReLU: GraphStar → Star
   S = rs.toStar();
   S_out = relu.reach(S, method, ...);

   % After ReLU: Star → GraphStar
   new_V = reshape(S_out.V, [numNodes, numFeatures, S_out.nVar + 1]);
   rs = GraphStar(new_V, S_out.C, S_out.d, ...
                  S_out.predicate_lb, S_out.predicate_ub);

How affineMap Works
^^^^^^^^^^^^^^^^^^^^

For affine layers (FullyConnectedLayer, Conv2DLayer), the layer applies its
weight matrix to each column of ``V`` (center + each generator):

.. code-block:: matlab

   % Star.affineMap simplified:
   function s2 = affineMap(obj, W, b)
       new_V = W * obj.V;           % Apply W to center and all generators
       new_V(:,1) = new_V(:,1) + b; % Add bias to center only
       s2 = Star(new_V, obj.C, obj.d, obj.predicate_lb, obj.predicate_ub);
   end

For specialized set types (ImageStar), the layer reshapes V to apply operations
in the native domain (e.g., 2D convolution on each generator image), then
reshapes back. The constraints ``C``, ``d`` are **always preserved unchanged**
through affine operations -- this is a key property that makes Star sets
efficient.

Integration Strategies
-----------------------

**Strategy 1: Conversion-based (recommended for most cases)**

Implement ``toStar()`` and a static ``fromStar()``. Your set flows through
affine layers as-is (via your own ``affineMap``), and converts to Star for
activation layers. This is how ImageStar, VolumeStar, and GraphStar work.

**Strategy 2: Custom network class**

If your set requires extra context (like GNN needs the adjacency matrix),
create a new network class (like ``GNN.m``) that overrides the reach loop.
This lets you pass additional arguments to specialized layers while still
using standard activation layers via conversion.

**Strategy 3: Extend existing layers**

Add cases to existing layer ``reach()`` methods to handle your set type
natively. Only recommended for simple extensions.

Example: How GraphStar Was Added
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Created ``GraphStar.m`` in ``engine/set/`` with V shaped as ``(N, F, numBasis+1)``
2. Implemented ``toStar()`` (flatten N*F) and reconstruction from Star
3. Created ``GCNLayer.m`` and ``GINELayer.m`` with graph-aware ``reach()``
4. Created ``GNN.m`` that routes graph layers vs activation layers
5. No changes needed to ``Star.m``, ``ReluLayer.m``, or ``NN.m``

This modular approach means new set types don't require modifying existing code.

Reference Files
----------------

- ``Star.m`` (~2050 lines) -- the fundamental set type
- ``ImageStar.m`` -- 2D image specialization with toStar/fromStar
- ``GraphStar.m`` -- graph node feature specialization
- ``GNN.m`` (lines 227-273) -- the conversion pattern in reach()
