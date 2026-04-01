# Probabilistic Verification Theory

## Motivation

For large neural networks -- particularly semantic segmentation networks (SSNs) with millions of parameters and high-dimensional, pixel-level outputs -- computing exact or even approximate reachable sets can be computationally intractable. The exponential complexity of ReLU branching and the sheer output dimensionality (e.g., $256 \times 512 \times 19$ logits for a Cityscapes model) render deterministic methods impractical.

NNV3 integrates a **probabilistic, data-driven verification** approach based on conformal inference that provides an alternative: instead of verifying properties for *all* possible perturbations, it verifies them with *high probability* under a given distribution. This approach is **architecture-agnostic** (scales with inference cost, not network structure) and provides formal probabilistic guarantees.

## Conformal Inference Framework

### Core Idea

Given a collection of i.i.d. scalar random variables $\mathbf{M} = \{R_1, R_2, \ldots, R_m\}$ sampled from $R \sim \mathcal{D}$ (called **nonconformity scores**), conformal inference constructs a prediction interval for a new unseen draw $R^{\text{unseen}}$ from the same distribution.

The key result: if we sort the calibration scores $R_1 < R_2 < \ldots < R_m$ and select rank $\ell := \lceil(m+1)(1-\epsilon)\rceil$, then:

$$\Pr[R^{\text{unseen}} \leq R_\ell] \geq 1 - \epsilon$$

where the probability is over the joint randomness of both calibration and test samples.

### Double-Step Probabilistic Guarantee

The coverage level $\delta = \Pr[R^{\text{unseen}} \leq R_\ell]$ is itself a random variable following a $\text{Beta}(\ell, m+1-\ell)$ distribution. By appropriately tuning $m$ and $\ell$, the variance of $\delta$ can be made extremely small (e.g., $3.123 \times 10^{-8}$ for $m = 8000$, $\ell = 7999$).

NNV3 uses the **double-step (PAC) guarantee** formulation using the regularized incomplete beta function:

$$\Pr\left[\Pr[R^{\text{unseen}} \leq R_\ell] > 1 - \epsilon\right] > 1 - \text{betacdf}_{1-\epsilon}(\ell, m+1-\ell)$$

This is denoted as the $\langle \epsilon, \ell, m \rangle$ guarantee, where:
- $\epsilon$ is the **miscoverage level** (e.g., 0.0001 for 99.99% coverage)
- $\delta_1 = 1 - \epsilon$ is the **coverage level**
- $\delta_2 = 1 - \text{betacdf}_{1-\epsilon}(\ell, m+1-\ell)$ is the **confidence** of the guarantee

A key advantage of the double-step formulation is its **flexibility**: the strict relationship among $m$, $\ell$, and $\epsilon$ can be relaxed, allowing practitioners to tune these parameters to achieve desired coverage-confidence tradeoffs.

### Example

With $m = 100{,}000$ calibration samples and rank $\ell = 99{,}999$:

$$\Pr\left[\Pr[R^{\text{unseen}} < R_{99{,}999}] > 0.9999\right] > 0.9995$$

This guarantee holds for *any* calibration dataset sampled from distribution $\mathcal{D}$, regardless of the specific distribution shape.

## Probabilistic Reachability for Neural Networks

### Problem Formulation

Given a neural network $f : \mathbb{R}^{n_0} \to \mathbb{R}^n$, an input set $\mathbf{I} \subset \mathbb{R}^{n_0}$, and a sampling distribution $x \overset{\mathcal{W}}{\sim} \mathbf{I}$, probabilistic reachability constructs a set $\mathbf{R}_f^\epsilon(\mathcal{W}; \ell, m)$ such that:

$$x \overset{\mathcal{W}}{\sim} \mathbf{I} \;\Rightarrow\; \Pr\left[\Pr[f(x) \in \mathbf{R}_f^\epsilon(\mathcal{W}; \ell, m)] \geq 1 - \epsilon\right] \geq 1 - \text{betacdf}_{1-\epsilon}(\ell, m+1-\ell)$$

### Nonconformity Score Design

For an output $\text{vec}(y) = [y(1), \ldots, y(n)] \in \mathbb{R}^n$, the nonconformity score is defined as:

$$R_i^{\text{calib}} = \max\left(\frac{|y_i(1) - c(1)|}{\tau_1}, \ldots, \frac{|y_i(n) - c(n)|}{\tau_n}\right)$$

where $c \in \mathbb{R}^n$ is the center (average of training outputs) and $\tau_k$ are normalization factors computed from the training dataset. This score design yields a **hyper-rectangular reachable set** centered at $c$ with half-widths $\tau_k R_\ell^{\text{calib}}$.

## Scaling to Semantic Segmentation Networks

### The High-Dimensionality Challenge

The naive hyper-rectangular approach becomes overly conservative in high-dimensional output spaces because:

1. **Shape limitation**: Hyper-rectangles are poor approximations of the true output distribution shape in high dimensions
2. **Distribution space**: The space of distributions $\mathcal{Y}$ that can represent the calibration data grows dramatically with dimension, making conformal inference inherently more conservative

### Surrogate-Based Approach

To address these challenges, the NeurIPS 2025 paper (Hashemi et al.) introduces a **surrogate-based reachability** technique:

1. **Train a small ReLU surrogate** $g : \mathbb{R}^{n_0} \to \mathbb{R}^n$ to approximate the original network $f$ on the input set $\mathbf{I}$
2. **Compute deterministic reachable set** of $g$ using NNV's Star-set reachability: $\mathbf{R}_g(\mathbf{I})$
3. **Compute prediction error** $q(x) = f(x) - g(x)$ and apply the naive CI technique to bound this error with an $\langle \epsilon, \ell, m \rangle$ guaranteed hyper-rectangle $\mathbf{R}_q^\epsilon$
4. **Combine via Minkowski sum**: $\mathbf{R}_f^\epsilon = \mathbf{R}_g(\mathbf{I}) \oplus \mathbf{R}_q^\epsilon$

This yields two key improvements:
- The reachable set is no longer constrained to a hyper-rectangle (the surrogate's Star set captures the output shape)
- The CI calibration is performed on *prediction errors* rather than raw outputs, which are much smaller in magnitude

### Dimensionality Reduction via Deflative PCA

For SSNs with extremely high-dimensional outputs (e.g., $720 \times 960 \times 12$ for CamVid), even training the surrogate is challenging. The approach uses a **two-stage training**:

1. **Stage 1 (PCA via deflation)**: Learn the top $N \ll n$ principal directions of the output point cloud using a learning-based deflation algorithm. This maps high-dimensional logits to a low-dimensional representation $v \in \mathbb{R}^N$ via $v = A^\top \text{vec}(y)$
2. **Stage 2 (ReLU network)**: Train a small ReLU network $g_2 : \mathbb{R}^r \to \mathbb{R}^N$ mapping perturbation coefficients to the low-dimensional space

The final surrogate is $g(x^{\text{adv}}) = A \cdot g_2(\lambda)$, which is scalable to SSNs.

### Pixel-Level Robustness Classification

Once the reachable set over SSN logits is constructed as a Star set, each pixel $(i,j)$ is classified by projecting onto its logit components:

- **Robust**: The correct class logit's lower bound exceeds all other classes' upper bounds
- **Non-robust**: A different class dominates for all perturbations
- **Unknown**: Logit intervals overlap between classes (due to over-approximation conservatism)

The **Robustness Value (RV)** = percentage of pixels classified as robust.

## Comparison with Deterministic Methods

| Aspect | Deterministic (Star) | Probabilistic (CP) |
|--------|---------------------|-------------------|
| Guarantee type | Worst-case (all inputs) | Statistical ($\langle \epsilon, \ell, m \rangle$ guarantee) |
| Completeness | Sound (and complete for exact) | PAC coverage guarantee |
| Scalability | Limited by network size/depth | Scales with inference cost |
| Requirements | LP solver, symbolic analysis | Sampling + surrogate training |
| Output dimensionality | Struggles with high-dim outputs | Handles SSNs via PCA + CI |
| Distribution assumption | None (worst-case) | Requires sampling distribution $\mathcal{W}$ |
| Best for | Small-medium networks | Large, intractable networks (SSNs, object detection) |

## Comparison with Randomized Smoothing

An advantage of conformal inference over randomized smoothing (Fischer et al., 2021; Anani et al., 2024) is that CI guarantees are defined **directly on the base model** $\text{SSN}(x')$, while randomized smoothing operates on a noise-injected approximation $\text{SSN}(x' + \nu)$. To make randomized smoothing apply to the base model, one would need $\nu = 0$ (i.e., $\sigma = 0$), which collapses the perturbation ball to a singleton.

On the other hand, randomized smoothing provides guarantees over the **worst-case** input $x' \in \mathbf{B}_r(x)$, whereas conformal inference assumes a **prior distribution** $x' \overset{\mathcal{W}}{\sim} \mathbf{B}_r(x)$.

## Limitations

- Requires a **prior sampling distribution** $\mathcal{W}$ (unlike worst-case deterministic methods)
- Calibration set must be **recomputed per test point** (main computational cost)
- Surrogate-based technique increases runtime (PCA + training + deterministic reachability)
- For very high perturbation dimensions, falls back to the naive technique for scalability

## See Also

- {doc}`/user-guide/conformal-prediction` -- Python setup and NNV usage guide
- {doc}`/examples/semantic-segmentation` -- Semantic segmentation verification examples
- {doc}`/examples/medical-imaging` -- Medical imaging verification with probabilistic methods
- N. Hashemi, S. Sasaki, I. Oguz, M. Ma, T.T. Johnson, "Scaling Data-Driven Probabilistic Robustness Analysis for Semantic Segmentation Neural Networks", NeurIPS 2025
