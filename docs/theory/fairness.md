# Fairness Verification

## Motivation

As neural networks are deployed for high-stakes decisions (loan approvals, credit scoring, hiring), ensuring algorithmic fairness becomes critical. Unlike statistical auditing (which tests fairness on sampled inputs), NNV's **FairNNV** module provides formal verification of fairness properties over continuous input regions.

## Fairness Definitions

### Counterfactual Fairness

A model is **counterfactually fair** if predictions remain unchanged when sensitive attributes (e.g., gender, race) are altered while all other attributes stay the same.

Formally, for a network $f$, input $x$, and sensitive attribute index $s$:

$$\forall x \in \mathcal{X}, \forall s' \in \mathcal{S}: \; f(x) = f(x[s \leftarrow s'])$$

NNV verifies this by constructing an input set where the sensitive attribute varies over its full range and checking whether the output reachable set changes.

### Individual Fairness

**Individual fairness** (Lipschitz fairness) requires that similar individuals receive similar predictions:

$$\|f(x) - f(x')\| \leq L \cdot \|x - x'\|$$

NNV verifies this by computing the output reachable set for a neighborhood of inputs and checking that the output range is within the required bound.

## Verified Fairness (VF) Score

The **VF score** is a quantitative, formally grounded measure of fairness introduced by FairNNV:

- Partitions the input space into regions
- For each region, verifies whether the fairness property holds
- The VF score is the fraction of regions for which fairness is verified

Unlike statistical metrics (which have sampling error), the VF score is **formally sound**: if a region is marked as fair, the property provably holds for ALL inputs in that region.

## Methodology

1. **Define sensitive attributes** and their ranges
2. **Construct input sets** where sensitive attributes vary
3. **Compute reachable outputs** using Star-based reachability
4. **Check fairness property** by comparing output sets across sensitive attribute values
5. **Compute VF score** over the input space

## Advantages over Statistical Auditing

| Approach | Guarantee | Coverage | False negatives |
|----------|-----------|----------|-----------------|
| Statistical auditing | None (sampling-based) | Tested samples only | Possible |
| FairNNV (reachability) | Formal proof | All inputs in region | None |

FairNNV guarantees hold for ALL inputs within the verified region, not just the samples tested.

## References

- A.M. Tumlin, D. Manzanas Lopez, P. Robinette, Y. Zhao, T. Derr, T.T. Johnson, "FairNNV: The Neural Network Verification Tool For Certifying Fairness," ICAIF 2024
