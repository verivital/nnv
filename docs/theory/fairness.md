# Fairness Verification (FairNNV)

## Motivation

As neural networks are increasingly deployed for high-stakes decisions -- loan approvals, credit scoring, hiring, fraud detection -- ensuring algorithmic fairness becomes critical. Empirical fairness metrics (disparate impact, equalized odds) operate over finite test sets and may fail to expose systematic bias. FairNNV provides **formal verification** of fairness properties over continuous input regions, offering guarantees that hold for *all* inputs within specified bounds, not just sampled ones.

FairNNV is integrated into NNV3 as a module that leverages the Star-set reachability engine to verify individual and counterfactual fairness, producing a quantitative **Verified Fairness (VF) score**.

## Fairness Definitions

### Counterfactual Fairness

A model is **counterfactually fair** if its prediction remains unchanged when the sensitive attribute (e.g., gender, race) is altered while all other features are fixed.

Formally, for a predictor $\hat{Y}$, non-sensitive attributes $X$, and sensitive attributes $A$:

$$\hat{Y}(X^{(i)}, a) = \hat{Y}(X^{(i)}, a') \quad \text{for any } a, a' \in A \text{ s.t. } a \neq a'$$

FairNNV verifies this by **flipping the sensitive attribute** to its complementary value (e.g., male $\to$ female) and checking whether the classification output changes. If the reachable output set under the flipped attribute still maps to the same class, the sample is certified counterfactually fair.

### Individual Fairness

A model is **individually fair** if similar individuals receive similar predictions, regardless of sensitive attributes. This is a stronger notion than counterfactual fairness because it also allows non-sensitive features to vary within a bounded neighborhood.

Formally, given a distance metric $d(\cdot, \cdot)$ and perturbation radius $\epsilon$:

$$\Phi_{\text{IF}}(x, x', y, y', \epsilon) \stackrel{\text{def}}{=} \left(\|X_x - X_{x'}\| \leq \epsilon \;\wedge\; A_{x'} \in \text{alter}(A_x)\right) \rightarrow (y = y')$$

where $X$ are non-sensitive attributes, $A$ are sensitive attributes, and $\text{alter}(A_x)$ denotes the set of values obtained by altering the sensitive attribute.

FairNNV verifies this by:
1. Flipping the sensitive attribute(s)
2. Applying bounded perturbations $\epsilon$ to non-sensitive numerical features
3. Computing the reachable output set via Star-set reachability
4. Checking if the output class remains invariant over the entire perturbed region

## Connection to Adversarial Robustness

FairNNV models fairness specifications **analogously to adversarial robustness**. In robustness verification, we ask: "does the classification change under input perturbation?" In fairness verification, we ask the same question but with a specific perturbation structure:

- **Counterfactual fairness**: Perturbation is restricted to the sensitive attribute dimension only ($\epsilon = 0$ for all other features)
- **Individual fairness**: Perturbation includes both the sensitive attribute flip and bounded variation ($\epsilon > 0$) in non-sensitive features

This means NNV's existing Star-set reachability infrastructure (exact and approximate methods) can be directly applied to fairness verification without architectural changes.

## Verified Fairness (VF) Score

The **Verified Fairness score** quantifies the proportion of inputs for which a model is certifiably fair:

$$\text{VF} = \frac{1}{n} \sum_{i=1}^{n} F_i \quad \text{where} \quad F_i = \begin{cases} 1 & \text{if CF} \wedge \text{IF} \to \text{fair} \\ 0 & \text{otherwise} \end{cases}$$

where $n$ is the number of samples in the verification set. The VF score is analogous to the Certifiable Robustness score used in adversarial robustness evaluation, providing a formally grounded measure rather than an empirical estimate.

## Verification Pipeline

The FairNNV pipeline operates as follows:

1. **Load model**: Import trained ONNX classifier into NNV via `matlab2nnv`
2. **Construct fairness perturbation**: For each test sample, create a Star set encoding the sensitive attribute flip (counterfactual) and/or bounded perturbation (individual fairness)
3. **Compute reachable set**: Propagate the Star set through the network using exact reachability analysis
4. **Check specification**: Verify that the reachable output set maps entirely to the same class as the original prediction
5. **Report**: Fair (verified), Not Fair (violation found), or Unknown (over-approximation inconclusive)

## Adversarial Debiasing Analysis

FairNNV also enables evaluation of **bias mitigation techniques**. The ICAIF 2024 paper compares original models with adversarially debiased models (using the AIF360 toolkit) and reveals a critical insight:

While adversarial debiasing **empirically improves** fairness metrics (disparate impact, equal opportunity difference, average odds difference), the **VF scores for debiased models are generally lower** than for original models. This discrepancy suggests that adversarial debiasing may address fairness on observed data but fail to generalize to inputs outside the training distribution.

This finding highlights the importance of integrating formal verification into fairness evaluation: **relying solely on empirical metrics may be misleading**, as they cannot capture bias in low-density regions of the input space that formal verification exhaustively explores.

## Scalability Considerations

FairNNV uses **exact reachability** via Star sets, providing sound and complete guarantees. However, each ReLU operation requires splitting the Star set, leading to worst-case $O(2^N)$ complexity for an $N$-neuron network. In practice:

- Counterfactual fairness verification is fast (< 1 second per sample) because the perturbation space is highly constrained (only the sensitive attribute changes)
- Individual fairness verification scales with $\epsilon$ and model size, ranging from seconds to minutes per sample
- For larger models, approximate reachability (`approx-star`) can be used as a tractable alternative, trading completeness for scalability

## Limitations and Future Directions

- Currently focuses on **single sensitive attribute** perturbation; multi-attribute fairness verification is planned
- Limited to **classification** models (regression fairness requires different similarity metrics)
- Only evaluates **individual and counterfactual fairness**; group fairness definitions (demographic parity, equalized odds) are not yet supported as verification specifications
- Adversarial debiasing analysis reveals a gap between empirical and formal fairness -- further investigation into why debiased models fail formal verification is needed

## See Also

- {doc}`/application-domains` -- Financial predictions domain
- {doc}`/user-guide/verification-methods` -- Verification methods used by FairNNV
- A.M. Tumlin, D. Manzanas Lopez, P. Robinette, Y. Zhao, T. Derr, T.T. Johnson, "FairNNV: The Neural Network Verification Tool For Certifying Fairness", ICAIF 2024. [DOI: 10.1145/3677052.3698677](https://doi.org/10.1145/3677052.3698677)
