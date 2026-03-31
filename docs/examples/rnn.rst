Recurrent Neural Networks
=========================

.. rst-class:: lead

   Verify robustness of RNN classifiers on sequential data with
   variable-length time horizons. This tutorial covers building RNNs
   from weight matrices and running sequence verification.

----

What You Will Learn
-------------------

- How to construct an RNN in NNV from weight matrices (input kernel, recurrent kernel, bias)
- How to build a hybrid architecture (RecurrentLayer + FullyConnectedLayers)
- How to verify classification robustness over time sequences
- How variable time horizons affect robustness guarantees

Background
----------

Recurrent neural networks process sequential data by maintaining a hidden
state that accumulates information over time. At each time step, the RNN
applies:

.. math::

   h_t = \sigma(W_i x_t + W_h h_{t-1} + b_h)

where :math:`W_i` is the input kernel, :math:`W_h` is the recurrent kernel,
and :math:`\sigma` is the activation function. After processing the full
sequence, a feedforward network classifies the final hidden state.

NNV verifies that adversarial perturbations to the input sequence cannot
change the classification at any time horizon.

Step 1: Construct the RNN from Weights
-----------------------------------------

.. code-block:: matlab

   % Load pre-trained weight matrices
   load('simple_rnn.mat');   % Contains: kernel, recurrent_kernel, bias
   load('dense.mat');        % Contains: weights and biases for FC layers

   % Build the RecurrentLayer
   rnn_params = struct;
   rnn_params.Wi = double(simple_rnn.kernel);            % Input kernel (input_dim x hidden_dim)
   rnn_params.Wh = double(simple_rnn.recurrent_kernel);  % Recurrent kernel (hidden_dim x hidden_dim)
   rnn_params.bh = double(simple_rnn.bias);              % Hidden bias (hidden_dim x 1)
   rnn_params.Wo = eye(size(rnn_params.Wh, 1));          % Output = hidden state (identity)
   rnn_params.fh = 'poslin';    % Hidden activation: ReLU
   rnn_params.fo = 'purelin';   % Output activation: linear

   L1 = RecurrentLayer(rnn_params);

   % Build feedforward classification layers after the RNN
   L2 = FullyConnectedLayer(dense.W1, dense.b1);
   L3 = ReluLayer();
   L4 = FullyConnectedLayer(dense.W2, dense.b2);
   L5 = ReluLayer();
   L6 = FullyConnectedLayer(dense.W3, dense.b3);
   % ... (continue for all dense layers)

   % Assemble the full network
   net = NN({L1, L2, L3, L4, L5, L6});

Step 2: Load Test Sequences
------------------------------

.. code-block:: matlab

   % Load test data
   load('points.mat');   % Contains sequence data

   % Select test sequences
   M = 5;   % Number of sequences to verify
   test_sequences = points(1:M, :)';   % Transpose to column vectors

Step 3: Verify Sequence Robustness
------------------------------------

For each test sequence, verify robustness at multiple time horizons.
A longer time horizon means more input steps can be perturbed, making
verification harder:

.. code-block:: matlab

   epsilon = 0.01;                    % Perturbation magnitude per time step
   time_horizons = [5, 10, 15, 20];  % Number of time steps to verify

   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';

   results = zeros(M, length(time_horizons));

   for k = 1:M
       x = test_sequences(:, k);

       for t = 1:length(time_horizons)
           Tmax = time_horizons(t);

           % Create input sequence: replicate the input point Tmax times
           % Each time step receives the same input, but perturbation
           % is applied independently at each step
           input_points = repmat(x, 1, Tmax);

           % Get true label by evaluating the unperturbed sequence
           y = net.evaluate(input_points);
           [~, target] = max(y);

           % Verify: can any ±epsilon perturbation at any time step
           % change the classification?
           result = net.verify_sequence_robustness(input_points, epsilon, target, reachOptions);

           results(k, t) = result;

           fprintf('Sequence %d, Tmax=%d: %s\n', k, Tmax, ...
               {'NOT ROBUST', 'ROBUST', 'UNKNOWN'}{result + 1});
       end
   end

Understanding the Results
--------------------------

.. code-block:: matlab

   % Display results table
   fprintf('\n         ');
   for t = time_horizons
       fprintf('T=%-4d ', t);
   end
   fprintf('\n');
   for k = 1:M
       fprintf('Seq %d:   ', k);
       for t = 1:length(time_horizons)
           labels = {'  X  ', '  OK ', '  ?  '};
           fprintf('%-6s ', labels{results(k,t) + 1});
       end
       fprintf('\n');
   end

Typical patterns:

- **Short horizons (T=5)**: Most sequences verified robust (small input space)
- **Medium horizons (T=10-15)**: Some sequences become unknown or not robust
- **Long horizons (T=20)**: Larger input space makes verification harder;
  more sequences may return unknown due to over-approximation accumulation

.. admonition:: Variable-Length Sequences (NNV3)

   NNV3 extends RNN verification beyond fixed-length sequences to handle
   variable-length inputs where T ranges from T_min to T_max. The reachable
   set is the union over all possible time horizons. This is useful for
   applications like NLP or sensor data where sequence length varies at runtime.

Source Files
^^^^^^^^^^^^

- `NN/RNN/ <https://github.com/verivital/nnv/tree/master/code/nnv/examples/NN/RNN>`_
