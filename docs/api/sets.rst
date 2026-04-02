Set Types
=========

.. rst-class:: lead

   NNV's set representations for reachability analysis. Each type captures
   a region of possible values with linear constraints on predicate variables.

----

Star
----

The fundamental polytopic set representation.

**Constructors:**

.. code-block:: matlab

   S = Star(lb, ub)                              % From box bounds
   S = Star(V, C, d)                             % From basis matrix and constraints
   S = Star(V, C, d, pred_lb, pred_ub)           % With predicate bounds

**Key Properties:** ``V`` (basis matrix), ``C``, ``d`` (constraints), ``dim``, ``nVar``,
``predicate_lb``, ``predicate_ub``, ``state_lb``, ``state_ub``

**Methods:**

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``S2 = S.affineMap(W, b)``
     - Apply affine transformation: y = Wx + b
   * - ``[lb, ub] = S.getRanges()``
     - Compute element-wise bounds via LP
   * - ``S2 = S.intersectHalfSpace(H, g)``
     - Intersect with half-space Hx <= g
   * - ``bool = S.contains(x)``
     - Check if point x is in the set
   * - ``bool = S.isEmptySet()``
     - Check if set is empty (infeasible constraints)
   * - ``S3 = S.MinkowskiSum(S2)``
     - Minkowski sum with another Star
   * - ``points = S.sample(N)``
     - Sample N random points from the set
   * - ``S.plot()``
     - Visualize the set (2D/3D projections)
   * - ``box = S.getBox()``
     - Get bounding box
   * - ``Z = S.getZono()``
     - Get zonotope outer approximation
   * - ``IS = S.toImageStar(h, w, c)``
     - Convert to ImageStar with given dimensions
   * - ``P = S.toPolyhedron()``
     - Convert to Polyhedron object

ImageStar
---------

Star set for 2D multi-channel images (H x W x C).

**Constructors:**

.. code-block:: matlab

   IS = ImageStar(IM, LB, UB)            % From image + perturbation bounds
   IS = ImageStar(lb_image, ub_image)     % From lower/upper bound images

**Key Methods:** Same as Star plus ``toStar()``, ``estimateRanges()``

**Key Properties:** ``height``, ``width``, ``numChannel``, ``im_lb``, ``im_ub``

VolumeStar
----------

Star set for 3D volumetric data (H x W x D x C).

**Constructors:**

.. code-block:: matlab

   VS = VolumeStar(VOL, LB, UB)          % From volume + perturbation bounds

**Key Methods:** Same interface as ImageStar, for 4D tensors.

**Key Properties:** ``height``, ``width``, ``depth``, ``numChannel``

GraphStar
---------

Star set for graph node features (N x F).

**Constructors:**

.. code-block:: matlab

   GS = GraphStar(NF, LB, UB)            % Node features + perturbation bounds

**Key Properties:** ``V`` (3D: numNodes x numFeatures x numBasis), ``numNodes``, ``numFeatures``

Zono (Zonotope)
---------------

Centrally symmetric polytope defined by center and generators.

**Constructor:**

.. code-block:: matlab

   Z = Zono(c, V)    % c: center, V: generator matrix

**Key Methods:** ``affineMap``, ``MinkowskiSum``, ``getBox``, ``toStar``

Box
---

Axis-aligned hyper-rectangle.

**Constructor:**

.. code-block:: matlab

   B = Box(lb, ub)

**Key Methods:** ``toStar()``, ``center``, ``generators``

HalfSpace
----------

Linear constraint Gx <= g, used for safety specifications.

**Constructor:**

.. code-block:: matlab

   HS = HalfSpace(G, g)

**Usage:** Passed to ``verify_specification()`` to define the unsafe output region.
