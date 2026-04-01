Fairness Verification (FairNNV)
===============================

.. rst-class:: lead

   Verify that neural network classifiers make fair decisions regardless
   of sensitive attributes like gender, race, or age. This tutorial covers
   counterfactual and individual fairness verification on financial
   benchmark datasets.

----

What You Will Learn
-------------------

- How to verify counterfactual fairness (predictions invariant under sensitive attribute flip)
- How to verify individual fairness (similar individuals get similar predictions)
- How to compute the Verified Fairness (VF) score
- How adversarial debiasing affects formal fairness guarantees

See :doc:`/theory/fairness` for the theoretical foundations.

Prerequisites
^^^^^^^^^^^^^

- Trained ONNX classifier (e.g., Adult Census income predictor)
- NNV installed with FairNNV module

Step 1: Load Model and Data
-----------------------------

.. code-block:: matlab

   % Load trained Adult Census classifier from ONNX
   dlnet = importNetworkFromONNX('adult_census_model.onnx', InputDataFormats='BC');
   net = matlab2nnv(dlnet);

   % Load test data (features + labels)
   load('adult_census_test.mat');  % X_test, y_test
   % Identify sensitive attribute index (e.g., gender = column 9)
   sensitive_idx = 9;

Step 2: Counterfactual Fairness Verification
----------------------------------------------

Flip the sensitive attribute while keeping all other features fixed:

.. code-block:: matlab

   n_samples = 100;
   cf_results = zeros(n_samples, 1);

   reachOptions = struct;
   reachOptions.reachMethod = 'approx-star';

   for i = 1:n_samples
       x = X_test(i, :)';
       target = y_test(i);

       % Create counterfactual: flip the sensitive attribute
       x_cf = x;
       x_cf(sensitive_idx) = 1 - x(sensitive_idx);  % binary flip

       % Create Star set spanning both original and counterfactual
       lb = min(x, x_cf);
       ub = max(x, x_cf);
       input_set = Star(lb, ub);

       % Verify: does the classification change?
       result = net.verify_robustness(input_set, reachOptions, target);
       cf_results(i) = result;  % 1=fair, 0=unfair, 2=unknown
   end

   cf_vf = sum(cf_results == 1) / n_samples * 100;
   fprintf('Counterfactual VF Score: %.1f%%\n', cf_vf);

Step 3: Individual Fairness Verification
------------------------------------------

Flip the sensitive attribute AND perturb non-sensitive features by epsilon:

.. code-block:: matlab

   epsilon_values = [0.02, 0.03, 0.05, 0.07, 0.10];

   for e = 1:length(epsilon_values)
       epsilon = epsilon_values(e);
       if_results = zeros(n_samples, 1);

       for i = 1:n_samples
           x = X_test(i, :)';
           target = y_test(i);

           % Perturb: flip sensitive + bounded perturbation on non-sensitive
           lb = x - epsilon;
           ub = x + epsilon;

           % Flip sensitive attribute
           lb(sensitive_idx) = min(x(sensitive_idx), 1 - x(sensitive_idx));
           ub(sensitive_idx) = max(x(sensitive_idx), 1 - x(sensitive_idx));

           input_set = Star(lb, ub);
           result = net.verify_robustness(input_set, reachOptions, target);
           if_results(i) = result;
       end

       if_vf = sum(if_results == 1) / n_samples * 100;
       fprintf('IF VF Score (eps=%.2f): %.1f%%\n', epsilon, if_vf);
   end

Experimental Results
---------------------

FairNNV was evaluated on three fairness benchmark datasets:

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 20 20

   * - Dataset
     - Model
     - Sensitive Attr.
     - CF VF (%)
     - IF VF at eps=0.02 (%)
   * - Adult Census
     - AC-1 (16-8)
     - Gender
     - 89
     - ~85
   * - Adult Census
     - AC-3 (50)
     - Gender
     - 87
     - ~65
   * - German Credit
     - GC-2
     - Gender
     - 77
     - ~80
   * - Bank Marketing
     - BM-1
     - Age
     - 89
     - ~90

Key observations:

- **Counterfactual fairness**: High VF scores (74--89%) with verification times
  under 0.03s per sample
- **Individual fairness**: VF degrades as epsilon increases; larger models show
  steeper decline and higher verification cost
- **Debiasing paradox**: Empirically debiased models (via AIF360 adversarial
  debiasing) show *lower* VF scores than originals, suggesting formal
  verification captures unfairness that statistical metrics miss

Adapting to Your Own Classifier
---------------------------------

1. Train your classifier on tabular data with identified sensitive attributes
2. Export to ONNX format
3. Identify the sensitive attribute column index
4. Choose epsilon values appropriate for your domain
5. Run counterfactual and individual fairness verification

Source Files
^^^^^^^^^^^^

- `FairNNV code <https://zenodo.org/records/13910015>`_
