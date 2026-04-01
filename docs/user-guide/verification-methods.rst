Verification Methods
====================

.. rst-class:: lead

   NNV provides multiple reachability methods that trade off between precision,
   completeness, and computational cost. This guide covers all methods,
   configuration options, and how to interpret results.

----

Core Workflow
-------------

Every NNV verification follows three steps:

1. **Define input set** -- the region of inputs to verify
2. **Compute reachable set** -- propagate the input set through the network
3. **Check property** -- test if the reachable set satisfies a safety specification

.. code-block:: matlab

   % Step 1: Input set
   input_set = Star(lb, ub);

   % Step 2: Reachability
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';
   output_sets = net.reach(input_set, reachOptions);

   % Step 3: Property check
   unsafe_region = HalfSpace(G, g);  % Gx <= g defines the unsafe region
   result = verify_specification(output_sets, unsafe_region);

Reachability Methods
--------------------

exact-star
^^^^^^^^^^

**Sound and complete** -- produces the exact reachable set by splitting
Star sets at each ReLU neuron where the input range crosses zero.

.. code-block:: matlab

   reachOptions.reachMethod = 'exact-star';

- **Output**: Multiple disjoint Star sets whose union is the exact reachable set
- **Precision**: Exact (no over-approximation)
- **Cost**: Exponential in the number of neurons that require splitting
- **Best for**: Small networks, critical safety-of-life applications, or when you
  need provably complete results (can prove both safety AND unsafety)

.. admonition:: Soundness

   If ``exact-star`` reports SAFE, the property provably holds for ALL inputs
   in the input set. If it reports UNSAFE, there provably EXISTS an input that
   violates the property.

approx-star
^^^^^^^^^^^^

**Sound but incomplete** -- over-approximates the reachable set using a
single Star set per layer, avoiding exponential splitting.

.. code-block:: matlab

   reachOptions.reachMethod = 'approx-star';

   % Optional: control precision with relaxFactor (0 to 1)
   reachOptions.relaxFactor = 0;    % Full LP optimization (most precise, slower)
   reachOptions.relaxFactor = 0.5;  % Partial optimization
   reachOptions.relaxFactor = 1;    % No LP optimization (fastest, least precise)

- **Output**: Single over-approximating Star set
- **Precision**: Controlled by ``relaxFactor`` (0 = most precise, 1 = fastest)
- **Cost**: Polynomial (no splitting)
- **Best for**: Medium to large networks where exact analysis is intractable
- **Default method** -- good balance of precision and speed

.. admonition:: Soundness

   If ``approx-star`` reports SAFE, the property provably holds. If it
   reports UNKNOWN, the result is inconclusive (the over-approximation may
   be too coarse).

approx-zono
^^^^^^^^^^^^

**Sound but incomplete** -- uses zonotope abstraction for fast over-approximation.

.. code-block:: matlab

   reachOptions.reachMethod = 'approx-zono';

- **Output**: Single zonotope (convertible to bounding box)
- **Precision**: Lower than approx-star, but faster
- **Cost**: Very low -- no LP solving required
- **Best for**: Quick initial analysis, very large networks, or when used as
  a pre-filter before more precise methods

abs-dom
^^^^^^^

**Sound but incomplete** -- abstract domain-based reachability.

.. code-block:: matlab

   reachOptions.reachMethod = 'abs-dom';

- Uses abstract interpretation techniques
- Comparable to zonotope methods in precision
- Useful for specific layer configurations

ODE Methods (for Neural ODE blocks)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For networks containing ``ODEblockLayer``, additional methods are available
to handle the continuous-time dynamics:

.. code-block:: matlab

   % Direct matrix exponential (for linear ODE blocks)
   reachOptions.reachMethod = 'direct';

   % Krylov subspace method (efficient for high-dimensional linear systems)
   reachOptions.reachMethod = 'krylov';

   % Numerical ODE integration (for nonlinear ODE blocks)
   reachOptions.reachMethod = 'ode45';

Method Selection Guide
----------------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 15 35

   * - Method
     - Sound
     - Complete
     - Speed
     - Recommended For
   * - ``exact-star``
     - Yes
     - Yes
     - Slow
     - Small networks, must prove/disprove
   * - ``approx-star``
     - Yes
     - No
     - Medium
     - General-purpose (default)
   * - ``approx-zono``
     - Yes
     - No
     - Fast
     - Quick screening, very large nets
   * - ``abs-dom``
     - Yes
     - No
     - Fast
     - Alternative to zonotope
   * - ``direct``
     - Yes
     - Yes
     - Fast
     - Linear Neural ODEs
   * - ``krylov``
     - Yes
     - Approx
     - Fast
     - High-dim linear ODEs
   * - ``ode45``
     - Yes
     - Approx
     - Medium
     - Nonlinear Neural ODEs

The reachOptions Structure
--------------------------

All verification options are passed through a MATLAB struct:

.. list-table::
   :header-rows: 1
   :widths: 25 15 15 45

   * - Field
     - Type
     - Default
     - Description
   * - ``reachMethod``
     - string
     - ``'approx-star'``
     - Reachability method (see above)
   * - ``relaxFactor``
     - double
     - ``0``
     - Relaxation for approx-star (0 = precise, 1 = fast)
   * - ``numCores``
     - int
     - ``1``
     - Number of parallel cores for computation
   * - ``lp_solver``
     - string
     - ``'linprog'``
     - LP solver: ``'linprog'``, ``'glpk'``, or ``'gurobi'``

**Example with all options:**

.. code-block:: matlab

   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';
   reachOptions.relaxFactor = 0;
   reachOptions.numCores = 4;       % Use 4 cores
   reachOptions.lp_solver = 'glpk'; % Use GLPK solver

   output_sets = net.reach(input_set, reachOptions);

Interpreting Results
--------------------

The ``verify_specification()`` function returns:

- **1** (SAFE / UNSAT): The property holds for ALL inputs in the input set.
  No output in the reachable set violates the safety specification.
- **0** (UNSAFE / SAT): A violation was found. There exists an input in the
  input set that produces an output violating the specification.
- **2** (UNKNOWN): The verification is inconclusive. The over-approximate
  reachable set intersects the unsafe region, but this may be due to
  approximation rather than a real violation.

.. code-block:: matlab

   result = verify_specification(output_sets, unsafe_region);
   switch result
       case 1
           disp('SAFE: Property verified for all inputs');
       case 0
           disp('UNSAFE: Counterexample found');
       case 2
           disp('UNKNOWN: Try a more precise method or smaller input set');
   end

Falsification (Pre-Verification)
---------------------------------

Before running full reachability analysis, NNV supports **random falsification**
-- sampling concrete inputs to quickly find counterexamples:

.. code-block:: matlab

   % Sample random points from the input set
   nSamples = 1000;
   for i = 1:nSamples
       x = input_set.sample(1);
       y = net.evaluate(x);
       if unsafe_region.contains(y)
           fprintf('Counterexample found at sample %d\n', i);
           break;
       end
   end

This is the strategy used in VNN-COMP submissions: try random sampling first,
then fall back to full reachability only if no counterexample is found.

Parallel Computation
--------------------

NNV supports multi-core parallel computation for reachability analysis:

.. code-block:: matlab

   reachOptions.numCores = 8;  % Use 8 parallel workers

This is most effective for:

- ``exact-star`` (parallelizes across Star set splits)
- Large networks where LP solving dominates computation time
- Multiple input sets processed simultaneously

Probabilistic Verification
---------------------------

For networks where full reachability is intractable (e.g., large semantic
segmentation networks), NNV offers **conformal prediction-based verification**:

.. code-block:: matlab

   reachOptions.reachMethod = 'approx-star';
   reachOptions.train_mode = 'Linear';   % or 'ReLU'
   reachOptions.coverage = 0.99;
   reachOptions.confidence = 0.99;

   result = verify_robustness_cp(net, input_set, reachOptions, target, numClasses);

See :doc:`conformal-prediction` for full setup and usage details, and
:doc:`/theory/probabilistic` for the theoretical foundations.
