ACAS Xu Collision Avoidance
===========================

.. rst-class:: lead

   Verify safety properties of the ACAS Xu airborne collision avoidance
   system -- the canonical benchmark for neural network verification.
   This tutorial walks through loading ONNX models, parsing VNNLIB
   specifications, and verifying formal safety properties.

----

What You Will Learn
-------------------

- How to load an ONNX neural network and convert it to NNV format
- How to parse VNNLIB specification files defining safety properties
- How to verify formal safety properties using ``verify_vnnlib()``
- How to interpret verification results and extract output bounds
- How to run batch verification across multiple networks

Background
----------

The Airborne Collision Avoidance System (ACAS Xu) is a family of 45 neural
networks designed to advise unmanned aircraft on collision avoidance maneuvers.
Each network maps 5 sensor inputs to 5 advisory outputs:

**Inputs:**

1. Distance to intruder (ft)
2. Angle to intruder relative to ownship heading (rad)
3. Angle to ownship relative to intruder heading (rad)
4. Speed of ownship (ft/s)
5. Speed of intruder (ft/s)

**Outputs** (advisory with lowest score is selected):

1. Clear-of-Conflict (COC)
2. Weak Left
3. Weak Right
4. Strong Left
5. Strong Right

Safety properties define constraints like "if the intruder is directly ahead,
the advisory should NOT be Clear-of-Conflict" -- critical for avoiding
mid-air collisions.

Prerequisites
-------------

- NNV installed and on the MATLAB path
- ACAS Xu ONNX networks in ``data/ACASXu/onnx/``
- VNNLIB property files in ``data/ACASXu/vnnlib/``

These files are included in the NNV repository.

Step 1: Load an ONNX Network
------------------------------

.. code-block:: matlab

   % Path to ACAS Xu networks
   onnx_dir = fullfile(nnvroot(), 'data', 'ACASXu', 'onnx');

   % Load one of the 45 networks (5x9 grid indexed by previous advisory and tau)
   net_file = fullfile(onnx_dir, 'ACASXU_run2a_1_1_batch_2000.onnx');

   % Import as MATLAB dlnetwork, then convert to NNV
   dlnet = importNetworkFromONNX(net_file, InputDataFormats='BC');
   net = matlab2nnv(dlnet);

   fprintf('Network loaded: %d layers, input size %s, output size %s\n', ...
       net.numLayers, mat2str(net.InputSize), mat2str(net.OutputSize));

Step 2: Load a VNNLIB Safety Property
---------------------------------------

VNNLIB is a standard format for specifying verification properties. Each file
defines input bounds and output constraints:

.. code-block:: matlab

   % Path to property files
   vnnlib_dir = fullfile(nnvroot(), 'data', 'ACASXu', 'vnnlib');

   % Load Property 3
   vnnlib_file = fullfile(vnnlib_dir, 'prop_3.vnnlib');
   [lb, ub, prop] = load_vnnlib(vnnlib_file);

   % lb, ub: input bounds (5x1 vectors)
   % prop: output constraints (cell array of HalfSpace objects)

   fprintf('Input bounds:\n');
   for i = 1:length(lb)
       fprintf('  x%d: [%.4f, %.4f]\n', i, lb(i), ub(i));
   end

.. admonition:: What's in a VNNLIB File?

   A VNNLIB file contains two parts:

   1. **Input constraints**: Bounds on each input variable (e.g., ``x1 >= 0.3``
      and ``x1 <= 0.5``). These define the region of the input space to verify.
   2. **Output constraints**: Linear inequalities on the outputs (e.g.,
      ``y1 <= y2`` meaning "the score for advisory 1 must not exceed advisory 2").
      These define the safety property.

   The ``load_vnnlib()`` function parses these into lower/upper bound vectors
   and HalfSpace constraint objects.

Step 3: Verify the Property
-----------------------------

.. code-block:: matlab

   % Create input set from the property bounds
   input_set = Star(lb, ub);

   % Configure verification
   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';

   % Method 1: Direct VNNLIB verification (simplest)
   result = net.verify_vnnlib(vnnlib_file, reachOptions);

   % Method 2: Manual verification (more control)
   output_sets = net.reach(input_set, reachOptions);
   for i = 1:length(prop)
       res = verify_specification(output_sets, prop{i});
       fprintf('Property clause %d: %s\n', i, ...
           {'UNSAFE (violated)', 'SAFE (holds)', 'UNKNOWN'}(res + 1));
   end

Step 4: Examine the Output Bounds
-----------------------------------

.. code-block:: matlab

   % Get the output reachable set bounds
   [out_lb, out_ub] = output_sets(1).getRanges();

   fprintf('\nOutput bounds:\n');
   advisories = {'COC', 'Weak Left', 'Weak Right', 'Strong Left', 'Strong Right'};
   for i = 1:5
       fprintf('  %12s: [%8.4f, %8.4f]\n', advisories{i}, out_lb(i), out_ub(i));
   end

   % The advisory with the LOWEST score is selected
   % If the COC (index 1) lower bound exceeds all other upper bounds,
   % then COC is never selected -- which is what Property 3 requires

Understanding Property 3
^^^^^^^^^^^^^^^^^^^^^^^^^^

Property 3 states: "If the intruder is directly ahead and moving towards the
ownship, the advisory should not be Clear-of-Conflict (COC)."

In VNNLIB terms:

- **Input region**: Intruder ahead (specific angle/distance bounds)
- **Safety requirement**: The COC output score should NOT be the minimum
  (i.e., COC should not be the selected advisory)

If verification returns **SAFE**, this means for ALL inputs within the
specified region, the network will never advise Clear-of-Conflict -- exactly
what we need for collision avoidance safety.

Step 5: Batch Verification
----------------------------

Verify all 45 networks against a property:

.. code-block:: matlab

   results = zeros(5, 9);     % 5x9 grid of ACAS Xu networks
   times = zeros(5, 9);

   for i = 1:5
       for j = 1:9
           % Load network
           filename = sprintf('ACASXU_run2a_%d_%d_batch_2000.onnx', i, j);
           filepath = fullfile(onnx_dir, filename);
           dlnet = importNetworkFromONNX(filepath, InputDataFormats='BC');
           net = matlab2nnv(dlnet);

           % Verify
           tic;
           results(i,j) = net.verify_vnnlib(vnnlib_file, reachOptions);
           times(i,j) = toc;

           fprintf('Network (%d,%d): %s (%.2fs)\n', i, j, ...
               {'UNSAFE','SAFE','UNKNOWN'}{results(i,j)+1}, times(i,j));
       end
   end

   fprintf('\nSummary: %d SAFE, %d UNSAFE, %d UNKNOWN\n', ...
       sum(results(:)==1), sum(results(:)==0), sum(results(:)==2));
   fprintf('Total time: %.1f seconds\n', sum(times(:)));

Interpreting Results
---------------------

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Result
     - Meaning
   * - **1 (SAFE)**
     - The property holds for ALL inputs in the specified region.
       No input can produce an unsafe advisory. This is a formal guarantee.
   * - **0 (UNSAFE)**
     - A violation was found. There exists an input in the specified region
       that produces an unsafe advisory.
   * - **2 (UNKNOWN)**
     - The approximate method cannot determine the result. The
       over-approximation overlaps with the unsafe region. Try ``exact-star``
       for a definitive answer (slower but complete).

Adapting to Your Own ONNX Model
---------------------------------

The same ONNX + VNNLIB workflow works for any network:

1. Export your network to ONNX format
2. Write VNNLIB specifications for your safety properties
3. Load and verify:

.. code-block:: matlab

   dlnet = importNetworkFromONNX('your_model.onnx', InputDataFormats='BC');
   net = matlab2nnv(dlnet);
   result = net.verify_vnnlib('your_property.vnnlib', reachOptions);

Source Files
^^^^^^^^^^^^

- `Tutorial/NN/ACAS Xu/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/Tutorial/NN/ACAS%20Xu>`_
- `NN/ACASXU/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/ACASXU>`_
