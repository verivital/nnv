# Star Set Reachability

## Definition

A **Star set** (generalized star) is a convex polytope defined by a center vector, a set of basis vectors (generators), and linear constraints on the predicate variables:

$$S = \{ x \in \mathbb{R}^n \mid x = c + V\alpha, \;\; C\alpha \leq d \}$$

where:
- $c \in \mathbb{R}^n$ is the center (anchor) point
- $V \in \mathbb{R}^{n \times m}$ is the basis matrix (columns are generators)
- $\alpha \in \mathbb{R}^m$ is the vector of predicate variables
- $C \in \mathbb{R}^{p \times m}$ and $d \in \mathbb{R}^p$ define the linear constraints

This representation is both expressive (can represent any bounded convex polytope) and efficient for affine operations (the constraints propagate unchanged through linear layers).

## Exact Reachability (exact-star)

The **exact-star** method computes the precise reachable set by handling each ReLU neuron that straddles zero by **splitting** the Star set into two:

For a single ReLU neuron $y = \max(0, x_i)$ where the input range of $x_i$ crosses zero:

**Case 1** ($x_i \geq 0$): $y = x_i$, add constraint $x_i \geq 0$

**Case 2** ($x_i < 0$): $y = 0$, add constraint $x_i < 0$ and set the corresponding generator column to zero

The exact reachable set is the **union** of all resulting Star sets:

$$R = S_1 \cup S_2 \cup \ldots \cup S_k$$

### Complexity

The number of output Star sets grows exponentially with the number of "crossing" neurons (those whose input range straddles zero). For a layer with $n$ crossing neurons, the worst case is $2^n$ output sets.

**Mitigation strategies:**
- Empty set detection: Prune infeasible branches early via LP
- Parallel computation: Process independent branches on multiple cores
- Input set partitioning: Smaller input sets yield fewer crossings

## Approximate Reachability (approx-star)

The **approx-star** method avoids exponential splitting by computing a single over-approximating Star set at each ReLU neuron.

For a crossing neuron with input range $[l, u]$ where $l < 0 < u$:

Instead of splitting, introduce a **new predicate variable** $\beta$ and constrain:

$$y \geq 0, \quad y \geq x_i, \quad y \leq \frac{u}{u - l}(x_i - l)$$

This creates a triangular over-approximation of the exact ReLU behavior. The resulting Star set has one additional predicate variable per crossing neuron.

### Relaxation Factor

The `relaxFactor` parameter (0 to 1) controls how much LP optimization is used:

- **relaxFactor = 0**: Full LP at every neuron to compute tight bounds (most precise, slowest)
- **relaxFactor = 1**: Use estimated bounds without LP (fastest, least precise)
- **0 < relaxFactor < 1**: LP is used for a fraction of neurons

## Affine Operations

Star sets propagate efficiently through affine (linear) layers. For a fully connected layer $y = Wx + b$:

$$S_{out} = \{ y \mid y = Wc + b + WV\alpha, \;\; C\alpha \leq d \}$$

The constraints ($C$, $d$) remain unchanged -- only the center and generators are transformed. This is a key efficiency advantage of the Star representation.

## Verification via Intersection

To check if a safety property $Gx \leq g$ (HalfSpace) is satisfied:

1. Compute the reachable output set(s) $R$
2. For each output Star set $S_i$, check if $S_i \cap \{x \mid Gx \leq g\}$ is empty
3. If all intersections are empty: **SAFE** (property holds)
4. If any intersection is non-empty with `exact-star`: **UNSAFE** (counterexample exists)
5. If any intersection is non-empty with `approx-star`: **UNKNOWN** (may be spurious)

The intersection check reduces to an LP feasibility problem, which is efficient to solve.

## References

- H.-D. Tran et al., "Star-Based Reachability Analysis for Deep Neural Networks," FM 2019
- H.-D. Tran et al., "NNV: A Tool for Verification of Deep Neural Networks and Learning-Enabled Autonomous CPS," CAV 2020
- D. Manzanas Lopez et al., "NNV 2.0: The Neural Network Verification Tool," CAV 2023
