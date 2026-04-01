# Probabilistic Verification Theory

## Motivation

For large neural networks (e.g., semantic segmentation with millions of parameters), computing exact or even approximate reachable sets can be computationally intractable. NNV's probabilistic verification module provides an alternative: **data-driven reachability** based on conformal inference that provides formal coverage guarantees without requiring full symbolic analysis.

## Conformal Inference Framework

Conformal prediction provides distribution-free prediction regions with guaranteed coverage. Given:

- A trained neural network $f$
- An input perturbation region $\mathcal{X}$
- A desired coverage level $1 - \epsilon$ (e.g., 0.99)
- A desired confidence level $1 - \delta$ (e.g., 0.99)

The algorithm:

1. **Sample** $N$ points from the input perturbation region
2. **Evaluate** the network on each sample
3. **Train a surrogate model** to approximate the network's behavior
4. **Compute nonconformity scores** measuring how well the surrogate fits
5. **Construct prediction regions** calibrated to achieve the desired coverage

### Coverage Guarantee

The resulting prediction region $\hat{R}$ satisfies:

$$\Pr\left[\Pr_{x \sim \mathcal{X}}[f(x) \in \hat{R}] \geq 1 - \epsilon\right] \geq 1 - \delta$$

This is a **PAC (Probably Approximately Correct)** guarantee: with probability at least $1 - \delta$, the prediction region covers at least $1 - \epsilon$ of the output distribution.

## Sample Size Computation

The required number of samples is determined by the `CP_specification()` function based on:

- $\epsilon$: miscoverage level (e.g., 0.01 for 99% coverage)
- $\delta$: failure probability (e.g., 0.01 for 99% confidence)
- $N_\text{dir}$: number of projection directions for dimensionality reduction
- $N$: calibration set size
- $N_s$: number of surrogate training samples

Larger values of $1 - \epsilon$ and $1 - \delta$ require more samples but provide stronger guarantees.

## Surrogate Models

NNV supports two surrogate architectures, trained in Python via PyTorch:

### Linear Surrogate

$$g(x) = Wx + b$$

- Fast training, few parameters
- Works well when the network is approximately linear in the perturbation region
- Trained by `Trainer_Linear.py`

### ReLU Surrogate

$$g(x) = W_2 \cdot \text{ReLU}(W_1 x + b_1) + b_2$$

- More expressive, captures nonlinear behavior
- Better for highly nonlinear networks
- Trained by `Trainer_ReLU.py`

The surrogate is trained to minimize the prediction error on sampled input-output pairs from the original network.

## Comparison with Deterministic Methods

| Aspect | Deterministic (Star) | Probabilistic (CP) |
|--------|---------------------|-------------------|
| Guarantee type | Worst-case (all inputs) | Statistical (most inputs) |
| Completeness | Sound (and complete for exact) | PAC coverage guarantee |
| Scalability | Limited by network size | Scales to large networks |
| Requirements | LP solver, symbolic analysis | Sampling + Python training |
| Best for | Small-medium networks | Large, intractable networks |

## References

- N. Hashemi, S. Sasaki, I. Oguz, M. Ma, T.T. Johnson, "Scaling Data-Driven Probabilistic Robustness Analysis for Semantic Segmentation Neural Networks," NeurIPS 2025
