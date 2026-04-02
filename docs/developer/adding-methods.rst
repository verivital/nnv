Adding New Reachability Methods
===============================

.. rst-class:: lead

   How to add a new reachability method to NNV's verification pipeline.

----

Where Methods Are Dispatched
----------------------------

Each layer's ``reach()`` method contains a switch on the method string:

.. code-block:: matlab

   % Example from a layer's reach() method
   function outputSet = reach(obj, inputSet, method, ...)
       switch method
           case 'exact-star'
               outputSet = obj.reach_exact(inputSet, ...);
           case 'approx-star'
               outputSet = obj.reach_approx(inputSet, ...);
           case 'approx-zono'
               outputSet = obj.reach_zono(inputSet, ...);
           case 'abs-dom'
               outputSet = obj.reach_absdom(inputSet, ...);
           case 'my-new-method'                          % Your addition
               outputSet = obj.reach_my_method(inputSet, ...);
       end
   end

Steps to Add a New Method
--------------------------

1. **Choose a method name** (e.g., ``'my-new-method'``)

2. **Implement in each relevant layer**: Add a case to the ``reach()`` switch
   in every layer type your method should support. At minimum: ``FullyConnectedLayer``,
   ``ReluLayer``, ``Conv2DLayer``.

3. **For affine layers**: Affine reachability is typically method-independent
   (the linear transformation is exact for any set type). You may reuse the
   existing affine implementation.

4. **For nonlinear layers (ReLU)**: This is where methods differ:

   - ``exact-star``: Splits the Star at each crossing neuron (exponential worst-case)
   - ``approx-star``: Introduces linear relaxation constraints (polynomial)
   - ``approx-zono``: Uses zonotope abstraction (no LP required)
   - Your method: implement your own over-approximation strategy

5. **Update validation** (optional): If your method requires new reachOptions
   fields, update ``NN.validate_reach_options()`` to recognize them.

How Existing Methods Differ
----------------------------

**exact-star (ReLU)**:

For each neuron where pre-activation bounds straddle zero, split into two
Star sets: one where the neuron is active (x >= 0, y = x) and one where
it's inactive (x < 0, y = 0). Returns multiple Star sets.

**approx-star (ReLU)**:

For each crossing neuron, add a single linear relaxation:
y >= 0, y >= x, y <= u/(u-l) * (x - l). Returns one Star set with
additional predicate variables. LP cost per neuron for tight bounds.

**approx-zono (ReLU)**:

Convert Star to zonotope, apply interval-based ReLU abstraction.
No LP solving required -- fastest but least precise.

**relax-star**:

Same as approx-star but only solves LPs for a fraction of crossing neurons
(controlled by relaxFactor). Remainder use interval bounds.
The ``relaxFactor`` parameter (0 to 1) controls the fraction that use
interval bounds instead of LP.
