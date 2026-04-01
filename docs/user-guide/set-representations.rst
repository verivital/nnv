Set Representations
===================

.. rst-class:: lead

   NNV uses set-based reachability analysis where input regions are represented
   as geometric sets that propagate through network layers. NNV provides 10
   set types for different data formats and precision requirements.

----

Overview
--------

.. list-table::
   :header-rows: 1
   :widths: 20 25 55

   * - Set Type
     - Data Format
     - When to Use
   * - ``Star``
     - Vectors (1D)
     - FFNNs, general-purpose verification
   * - ``ImageStar``
     - Images (H x W x C)
     - CNN verification, image robustness
   * - ``VolumeStar``
     - Volumes (H x W x D x C)
     - 3D medical images, video classification
   * - ``GraphStar``
     - Graph nodes (N x F)
     - GNN verification, power systems
   * - ``Zono``
     - Vectors (1D)
     - Fast approximate analysis, nonlinear systems
   * - ``ImageZono``
     - Images (H x W x C)
     - Fast approximate CNN analysis
   * - ``Box``
     - Vectors (1D)
     - Simple bound specifications, input regions
   * - ``HalfSpace``
     - Constraints
     - Safety property specifications
   * - ``SetTree``
     - Hierarchical
     - Temporal tracking in control systems
   * - ``Conversion``
     - Utility
     - Converting between set representations

Star
----

The **Star set** is NNV's fundamental representation -- a convex polytope defined
by a center, basis vectors, and linear constraints on the predicate variables.

**Mathematical Definition:**

.. math::

   S = \{ x \in \mathbb{R}^n \mid x = c + \sum_{i=1}^{m} \alpha_i v_i, \;\;
   C \alpha \leq d \}

where :math:`c` is the center, :math:`v_i` are basis vectors (generators),
:math:`\alpha \in \mathbb{R}^m` are predicate variables, and
:math:`C\alpha \leq d` are linear constraints.

**Constructors:**

.. code-block:: matlab

   % From bounds (creates a box-shaped Star)
   S = Star(lb, ub);

   % From basis matrix and constraints
   % V = [c, v1, v2, ..., vm]  (center + generators)
   S = Star(V, C, d);

   % With known predicate bounds
   S = Star(V, C, d, predicate_lb, predicate_ub);

**Key Methods:**

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``affineMap(W, b)``
     - Apply affine transformation: y = Wx + b
   * - ``getRanges()``
     - Compute element-wise [lb, ub] bounds via LP
   * - ``getVertices()``
     - Enumerate vertices of the polytope
   * - ``isEmptySet()``
     - Check if the set is empty (infeasible constraints)
   * - ``contains(x)``
     - Check if a point is contained in the set
   * - ``intersectHalfSpace(H, g)``
     - Intersect with half-space Hx <= g
   * - ``convexHull(S2)``
     - Compute convex hull with another Star
   * - ``MinkowskiSum(S2)``
     - Minkowski sum with another Star
   * - ``sample(N)``
     - Sample N random points from the set
   * - ``plot()``
     - Visualize the set (2D or 3D projections)

**Example:**

.. code-block:: matlab

   % Create a 2D Star set (diamond shape)
   lb = [-1; -1];
   ub = [1; 1];
   S = Star(lb, ub);

   % Apply an affine map
   W = [2 0; 0 3];
   b = [1; 1];
   S_mapped = S.affineMap(W, b);

   % Get output bounds
   [lb_out, ub_out] = S_mapped.getRanges();

ImageStar
---------

**ImageStar** extends the Star set concept to multi-channel 2D images,
enabling efficient verification of CNNs.

**Mathematical Definition:**

.. math::

   \Theta = \{ x \in \mathbb{R}^{H \times W \times C} \mid
   x = c + \sum_{i=1}^{m} \alpha_i v_i, \;\; C\alpha \leq d \}

where :math:`c` is the center image and :math:`v_i` are generator images
of the same dimensions.

**Constructors:**

.. code-block:: matlab

   % From an image and perturbation bounds
   IS = ImageStar(IM, LB, UB);
   % IM: center image (H x W x C)
   % LB: lower bound perturbation (H x W x C)
   % UB: upper bound perturbation (H x W x C)

   % L-inf perturbation around an image
   epsilon = 0.01;
   LB = -epsilon * ones(size(IM));
   UB = epsilon * ones(size(IM));
   IS = ImageStar(IM, LB, UB);

**Key Methods:**

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``affineMap(W, b)``
     - Apply affine transformation
   * - ``getRanges()``
     - Pixel-wise [lb, ub] bounds
   * - ``contains(image)``
     - Check if an image is in the set
   * - ``toStar()``
     - Convert to flattened Star set
   * - ``sample(N)``
     - Sample N random images from the set
   * - ``estimateRanges()``
     - Fast bound estimation (less precise than LP)

VolumeStar
----------

**VolumeStar** extends ImageStar to 4D data -- 3D spatial volumes (e.g., MRI,
CT scans) or video frames (height x width x channels x frames).

**Mathematical Definition:**

.. math::

   V = \{ x \in \mathbb{R}^{H \times W \times C \times F} \mid
   x = c + \sum_{i=1}^{m} \alpha_i v_i, \;\; P\alpha \leq q \}

where :math:`c` is the center volume and :math:`v_i` are generator volumes.

**Constructors:**

.. code-block:: matlab

   % From a volume and perturbation bounds
   VS = VolumeStar(VOL, LB, UB);
   % VOL: center volume (H x W x C x F)
   % LB: lower bound perturbation
   % UB: upper bound perturbation

**Use Cases:**

- Robustness verification of video classification networks (C3D, I3D)
- Verification of 3D medical image classifiers (MRI, CT)
- Handling both spatial AND temporal perturbations

GraphStar
---------

**GraphStar** represents perturbation regions over node features in
graph-structured data, enabling GNN verification.

**Constructors:**

.. code-block:: matlab

   % Node feature perturbation
   GS = GraphStar(NF, LB, UB);
   % NF: nominal node features (N x F matrix)
   % LB: lower bound perturbation per node/feature
   % UB: upper bound perturbation per node/feature

**Use Cases:**

- Verifying GCN/GINE predictions under input uncertainty
- Power system safety verification (IEEE bus benchmarks)
- Node feature robustness analysis

Zono (Zonotope)
---------------

A **Zonotope** is a centrally symmetric polytope defined by a center and
generators. Zonotopes are more efficient than Stars for certain operations
but less precise for ReLU networks.

**Mathematical Definition:**

.. math::

   Z = \{ x \mid x = c + \sum_{i=1}^{m} \alpha_i g_i, \;\;
   -1 \leq \alpha_i \leq 1 \}

**Constructors:**

.. code-block:: matlab

   Z = Zono(c, V);
   % c: center vector
   % V: generator matrix (columns are generators)

**Key Methods:**

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Method
     - Description
   * - ``affineMap(W, b)``
     - Apply affine transformation
   * - ``MinkowskiSum(Z2)``
     - Minkowski sum with another zonotope
   * - ``convexHull(Z2)``
     - Convex hull with another zonotope
   * - ``getBox()``
     - Over-approximate with a bounding box
   * - ``toStar()``
     - Convert to Star set representation

ImageZono
---------

**ImageZono** is the zonotope analogue of ImageStar for image data.
It is used with the ``approx-zono`` reachability method.

.. code-block:: matlab

   IZ = ImageZono(IM, LB, UB);

Box
---

A **Box** (axis-aligned hyper-rectangle) is the simplest set representation,
defined by lower and upper bounds per dimension.

.. code-block:: matlab

   B = Box(lb, ub);

   % Convert to Star for verification
   S = B.toStar();

   % Get center and generators
   c = B.center;
   G = B.generators;  % diagonal matrix

Boxes are commonly used for specifying input regions and for fast
over-approximation of more complex sets.

HalfSpace
----------

A **HalfSpace** defines a linear constraint :math:`Gx \leq g`, used
primarily for specifying safety properties:

.. code-block:: matlab

   % Safety property: output y1 <= 0
   HS = HalfSpace([1 0], 0);

   % Check if reachable set intersects the unsafe region
   result = verify_specification(output_sets, HS);

Set Conversions
---------------

NNV supports conversions between set types:

.. code-block:: matlab

   % Box to Star
   S = Star(lb, ub);       % Direct construction
   S = B.toStar();         % From Box object

   % Zonotope to Star
   S = Z.toStar();

   % Star to Box (over-approximation)
   [lb, ub] = S.getRanges();

   % ImageStar to Star (flattening)
   S = IS.toStar();

Choosing the Right Set Type
---------------------------

- **FFNNs with few inputs**: Use ``Star`` with ``exact-star`` for provably complete results
- **CNNs / image inputs**: Use ``ImageStar`` -- it preserves spatial structure through conv layers
- **3D / video inputs**: Use ``VolumeStar`` for Conv3D-based networks
- **GNNs**: Use ``GraphStar`` for node feature perturbations
- **Fast approximation**: Use ``Zono`` / ``ImageZono`` with ``approx-zono``
- **Simple bounds only**: Use ``Box`` for quick input specification
- **Safety properties**: Use ``HalfSpace`` to define the unsafe region
