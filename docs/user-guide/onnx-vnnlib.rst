ONNX & VNNLIB Support
=====================

.. rst-class:: lead

   NNV supports standard neural network verification formats: ONNX for
   model exchange and VNNLIB for property specifications, enabling
   participation in VNN-COMP benchmarks.

----

Loading ONNX Models
-------------------

NNV provides two approaches for loading ONNX models:

**Using MATLAB's importer + matlab2nnv (recommended):**

.. code-block:: matlab

   % Import ONNX as MATLAB dlnetwork
   dlnet = importNetworkFromONNX('model.onnx', InputDataFormats='BC');

   % Convert to NNV format
   net = matlab2nnv(dlnet);

**Using the dedicated onnx2nnv converter:**

.. code-block:: matlab

   net = onnx2nnv('model.onnx');

The ``InputDataFormats`` parameter specifies dimension ordering:

- ``'BC'`` -- Batch x Channels (for 1D feature inputs)
- ``'BCSS'`` -- Batch x Channels x Spatial x Spatial (for 2D image inputs)
- ``'BCSSS'`` -- Batch x Channels x Spatial x Spatial x Spatial (for 3D volumes)

Converting MATLAB Networks
--------------------------

The ``matlab2nnv()`` function converts MATLAB deep learning formats to NNV:

.. code-block:: matlab

   % From SeriesNetwork (sequential model)
   net = matlab2nnv(series_network);

   % From DAGNetwork (directed acyclic graph)
   net = matlab2nnv(dag_network);

   % From dlnetwork (modern MATLAB format)
   net = matlab2nnv(dl_network);

Loading from .mat files
-----------------------

.. code-block:: matlab

   net = load_NN_from_mat('network.mat');

Parsing VNNLIB Specifications
-----------------------------

VNNLIB is the standard format for specifying verification properties, used in
the VNN-COMP competition. NNV's ``load_vnnlib()`` parses these files:

.. code-block:: matlab

   [lb, ub, prop] = load_vnnlib('property.vnnlib');
   % lb: input lower bounds vector
   % ub: input upper bounds vector
   % prop: output property constraints (HalfSpace representations)

VNNLIB files specify:

- **Input constraints**: Bounds on each input dimension
- **Output properties**: Linear constraints on outputs (AND/OR combinations)

Full Verification Workflow
--------------------------

The typical ONNX + VNNLIB workflow:

.. code-block:: matlab

   % 1. Load network
   dlnet = importNetworkFromONNX('model.onnx', InputDataFormats='BCSS');
   net = matlab2nnv(dlnet);

   % 2. Load property
   [lb, ub, prop] = load_vnnlib('property.vnnlib');

   % 3. Create input set
   input_set = Star(lb, ub);

   % 4. Configure verification
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';

   % 5. Verify
   result = net.verify_vnnlib('property.vnnlib', reachOptions);

Alternatively, using the manual approach for more control:

.. code-block:: matlab

   % Compute reachable set
   output_sets = net.reach(input_set, reachOptions);

   % Check each output property
   for i = 1:length(prop)
       res = verify_specification(output_sets, prop{i});
       fprintf('Property %d: %s\n', i, ...
           {'UNSAFE', 'SAFE', 'UNKNOWN'}{res + 1});
   end

Exporting to VNNLIB
--------------------

NNV can also export verification specifications to VNNLIB format:

.. code-block:: matlab

   export2vnnlib(lb, ub, prop, 'output_spec.vnnlib');

VNN-COMP Integration
--------------------

NNV participates in the annual `VNN-COMP <https://sites.google.com/view/vnn2024>`_
competition. The competition runner script follows this pattern:

1. **Random falsification** (100 samples) to quickly find counterexamples
2. **Reachability analysis** if no counterexample found
3. **Multiple solver fallback** (try Gurobi, then GLPK, then linprog)

See the competition submissions in
`examples/Submission/VNN_COMP2024/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Submission/VNN_COMP2024>`_
for the full implementation.
