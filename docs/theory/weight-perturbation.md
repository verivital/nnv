# Weight Perturbation Verification (ModelStar)

## Motivation

Neural networks deployed in practice are subject to parameter uncertainty arising from quantization, numerical imprecision, and hardware faults. Computing-in-Memory (CiM) DNN accelerators, for example, program weights onto analog devices (e.g., memristors) where small device variations can cause significant performance drops. Bit-flip attacks can modify specific weights to alter predictions for targeted inputs while preserving overall accuracy. Even standard quantization from floating-point to fixed-point can introduce perturbations that accumulate across layers.

Unlike input perturbation verification (which asks "is the network robust to noisy inputs?"), weight perturbation verification asks: **"is the network robust to noisy parameters?"** This is critical for deploying neural networks on resource-constrained or unreliable hardware.

NNV3 introduces **ModelStar**, a star-set-based representation for interval-bounded weight perturbations that enables reachability analysis under simultaneous input and parameter uncertainty.

## Perturbation Model

### Linear Layers

ModelStar focuses on perturbations in the parameters of **linear layers** (fully-connected and convolutional). A linear layer computes:

$$L_i(x) = W_i \odot x + b_i$$

where $\odot$ is a linear operation (matrix product for FC layers, convolution for conv layers), $W_i$ is the weight tensor, and $b_i$ is the bias tensor.

### Perturbation Maps

Each perturbation $\alpha_j$ affects specific parameters in a layer. A **perturbation map** $M_{W_j}$ is a tensor with the same dimensions as $W_i$, where each entry is non-zero if the corresponding weight is affected by $\alpha_j$, and zero otherwise. Similarly, a bias perturbation map $M_{b_j}$ identifies affected biases.

The perturbed weight matrix is:

$$W' = W + \sum_{j=1}^{m} \alpha_j M_{W_j}$$

where each $\alpha_j$ is bounded: $\alpha_j^l \leq \alpha_j \leq \alpha_j^u$.

This formulation is flexible: a single perturbation can affect multiple weights (e.g., all weights in a row), and different perturbations can affect different subsets of weights. This generality enables modeling of diverse hardware error patterns.

ModelStar is a theoretical framework for star-set encoding of **any** perturbations in a parameterized linear layer, expressed as intervals (upper and lower bounds) on weights. This distinguishes it from approaches like interval neural networks (INNs), which use interval arithmetic to reduce the number of neurons and compute over-approximations of a network's output range. While INN-based methods have some overlap with ModelStar, they primarily serve as a verification method for input perturbations rather than a general parameter uncertainty framework. However, if the goal is to check whether the output sets of an INN remain safe given an input set, the INN's weight intervals can be used directly as weight perturbation specifications in ModelStar for verification.

### Connection to Quantization

A common source of weight perturbation is the **quantized compression** of models, where weights are stored at reduced precision and computation proceeds in floating-point after loading from memory. Such errors are naturally modeled as bounded perturbations on the weights, where each parameter lies within a specified interval determined by the quantization resolution.

For an $n$-bit quantization scheme, the maximum relative error induced by rounding (half of the least significant bit) corresponds to approximately:

| Bit Width | Max Relative Error |
|-----------|-------------------|
| 10-bit | 0.0488% |
| 9-bit | 0.0977% |
| 8-bit | 0.1953% |
| 7-bit | 0.3906% |

These values closely align with the perturbation magnitudes (0.05%, 0.1%, 0.2%, 0.4%) used in ModelStar's experimental evaluation, meaning the existing results already capture robustness against compression-induced quantization errors.

Note that if the concern is the propagation of quantization errors under **fixed-point or integer arithmetic** during inference (rather than compression-induced weight rounding), this requires modeling nonlinearities introduced by discrete arithmetic operations, which are not currently supported within ModelStar and remain part of planned future work.

## Reachability Analysis

### Single-Layer Perturbation (Exact)

When the input to a perturbed linear layer is a **singleton** (a set containing a single element -- i.e., no prior perturbations), the output reachable set is represented **exactly** by a Star set.

Given input $\{x\}$ and $m$ perturbations $\{\alpha_1, \ldots, \alpha_m\}$ with perturbation maps $\{M_{W_1}, \ldots, M_{W_m}\}$ and $\{M_{b_1}, \ldots, M_{b_m}\}$:

- **Center**: $c = W_i \odot x + b_i$ (the nominal output)
- **Basis tensors**: $v_j = M_{W_j} \odot x + M_{b_j}$ for each perturbation $\alpha_j$
- **Constraints**: $\alpha^l \leq \alpha \leq \alpha^u$ (box bounds on perturbation variables)

The resulting Star set $\langle c, V, P \rangle$ with $m$ generators exactly captures all possible outputs of the perturbed layer. No over-approximation is introduced at this stage.

### Multi-Layer Perturbation (Sound Over-Approximation)

When the input to a perturbed linear layer is already a Star set (from prior perturbations or input uncertainty), the new perturbation variables multiply with existing predicate variables, creating **nonlinear** terms. ModelStar handles this via **over-approximation** to preserve the linear Star set structure.

Given:
- Input Star set with $m$ existing generators
- Layer with $n$ neurons and $q$ new perturbations

The output reachable set is over-approximated by a Star set with $m + n + q$ generators. The additional $n$ generators arise from bounding the cross-terms between existing and new perturbation variables.

### Propagation Through the Network

The full ModelStar workflow:

1. **Before the perturbed layer**: Standard Star set propagation (affine layers are exact, ReLU uses approx-star)
2. **At the perturbed layer**: Apply the perturbation map formulation (exact for first perturbed layer, over-approximate for subsequent)
3. **After the perturbed layer**: Continue standard Star set propagation through remaining layers
4. **Safety check**: Intersect the final reachable set with the unsafe region

The reachable set starts as a singleton (single input image), remains a singleton through unperturbed layers, **expands** at the perturbed layer, and then propagates through subsequent layers with the expanded uncertainty.

## Safety Verification

For classification tasks, the unsafe region is the union of half-planes where an incorrect class score exceeds the correct class score:

$$\text{Unsafe} = \bigcup_{j \neq c} \{ y \mid y_j \geq y_c \}$$

where $c$ is the correct class. If the reachable set does not intersect the unsafe region, the network is **verified robust** against the specified weight perturbations. If it does intersect (under over-approximation), the result is **unknown**.

## Computational Complexity

The complexity of ModelStar reachability is governed by the number of generators in the Star set, which determines the cost of LP-based ReLU over-approximation in subsequent layers:

- **First perturbed layer**: Produces $m$ generators (one per perturbation variable). With one perturbation per weight, $m$ equals the number of weights in the layer.
- **Subsequent perturbed layers**: Each adds $n + q$ generators (where $n$ is the number of neurons and $q$ is the number of new perturbations), due to the over-approximation of cross-terms.
- **ReLU layers**: Each crossing neuron adds constraints via the approx-star relaxation, with cost scaling with the number of generators.

For layers with large parameter counts (e.g., fc_1 with 802,816 weights), execution time increases significantly. This motivates focusing verification on later, smaller layers or using selective perturbation maps that target specific weight subsets rather than perturbing all weights simultaneously.

## Current Limitations and Future Directions

- ModelStar currently supports **fully-connected and 2D convolutional layers**. Extension to other linear layers (batch normalization, 1D/3D convolution) is planned.
- Verification is performed **one input at a time** -- batch verification across multiple inputs is not yet supported.
- **Fixed-point/integer arithmetic** effects during inference (as opposed to weight compression errors) require modeling discrete arithmetic nonlinearities, which are not currently supported and remain future work.
- Multi-layer simultaneous perturbation is supported but yields lower verification rates due to accumulated over-approximation.

## See Also

- {doc}`/user-guide/verification-methods` -- Verification methods including weight perturbation
- {doc}`/user-guide/set-representations` -- Star set and ModelStar definitions
- {doc}`/application-domains` -- Industrial prognostics domain (Collins RUL with weight perturbation)
- M.U. Zubair, T.T. Johnson, K. Basu, W. Abbas, "Verification of Neural Network Robustness Against Weight Perturbations Using Star Sets", IEEE CAI 2025. [DOI: 10.1109/CAI64502.2025.00117](https://doi.org/10.1109/CAI64502.2025.00117)
