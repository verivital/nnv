# Graph Neural Network Reachability

## Motivation

Graph neural networks (GNNs) are increasingly used as fast, topology-aware surrogates in safety-critical domains such as power systems, molecular modeling, and traffic forecasting. Unlike feedforward networks that operate on flat vectors, GNNs compute over graph domains using **message-passing layers**, where nodes and edges carry feature vectors that are exchanged and transformed. Formally verifying GNNs requires lifting classical vector-based reachability analysis to graph-structured inputs where perturbations propagate through neighbor aggregation.

NNV's **GNNV** module addresses this by introducing **GraphStar sets** -- a generalization of Star sets that captures uncertainty over both node and edge features -- and implementing reachability analysis for two GNN architectures: GCN and GINE.

## GNN Architectures

### Graph Convolutional Network (GCN)

GCN layers aggregate and normalize features from neighboring nodes using a symmetric normalization of the adjacency matrix. Self-loops are added so each node also contributes its own representation:

$$H^{(\ell+1)} = \sigma\left(\tilde{D}^{-\frac{1}{2}} \tilde{A} \tilde{D}^{-\frac{1}{2}} H^{(\ell)} W^{(\ell)}\right)$$

where $\tilde{A} = A + I_N$ adds self-loops, $\tilde{D}$ is the corresponding degree matrix, $W^{(\ell)}$ is a trainable weight matrix, and $\sigma$ is an element-wise activation (ReLU).

### Graph Isomorphism Network with Edge Features (GINE)

GINE layers incorporate both node **and** edge features into each message-passing step, using a two-layer MLP with an internal ReLU:

**Aggregation** (per-edge messages summed at target nodes):

$$Z^{(\ell)} = R^\top\left(S H^{(\ell)} + E W_e^{(\ell)} + b_e^{(\ell)}\right)$$

**Combination** (two-layer MLP):

$$H^{(\ell+1)} = \sigma\left(\sigma\left(Z^{(\ell)} W_1^{(\ell)} + b_1^{(\ell)}\right) W_2^{(\ell)} + b_2^{(\ell)}\right)$$

where $S, R \in \{0,1\}^{M \times N}$ are the source and target incidence matrices, $E$ is the edge-feature matrix, and $W_e$ projects edge features into the node-feature dimension. The outer ReLU is omitted at the final layer.

The key distinction from GCN is that GINE explicitly uses edge information (e.g., line impedances in power grids), enabling verification under joint node and edge perturbations.

## GraphStar Sets

### Definition

GraphStar sets extend the Star set domain to graph-structured feature matrices. Two complementary abstractions capture uncertainty:

**NodeGraphStar** $\Theta_X = \langle C_X, B_X, P_X \rangle$ represents a set of node-feature matrices:

$$[\![\Theta_X]\!] = \left\{ X \;\middle|\; X = C_X + \sum_{i=1}^{m} \alpha_i B_i, \;\; \Gamma_X \alpha \leq \gamma_X \right\}$$

where $C_X \in \mathbb{R}^{N \times d_f}$ is the center, $\{B_1, \ldots, B_m\}$ are generator matrices in $\mathbb{R}^{N \times d_f}$, and $\Gamma_X \alpha \leq \gamma_X$ constrains the predicate variables.

**EdgeGraphStar** $\Theta_E = \langle C_E, B_E, P_E \rangle$ represents a set of edge-feature matrices:

$$[\![\Theta_E]\!] = \left\{ E \;\middle|\; E = C_E + \sum_{j=1}^{r} \beta_j B_j, \;\; \Gamma_E \beta \leq \gamma_E \right\}$$

Joint uncertainty is denoted $\Theta_\mathcal{G} = (\Theta_X, \Theta_E)$.

### Why GraphStar?

Unlike standard Star sets which operate on flat vectors, GraphStar sets maintain the **matrix structure** $\mathbb{R}^{N \times d_f}$ (or $\mathbb{R}^{M \times d_e}$) so that graph operations -- neighbor aggregation via $\tilde{A}$, source gathering via $S$, and target scattering via $R^\top$ -- can be applied directly to the center and generator matrices as standard matrix multiplications. This avoids flattening and reshaping, preserving the row-wise correspondence between features and graph entities.

### Construction from Perturbation Specs

For an $\ell_\infty$ node perturbation of magnitude $\epsilon_{\text{node}}$, the NodeGraphStar has center $C_X = X$ with one generator matrix per perturbed feature column, where each $B_i$ encodes the perturbation range across all $N$ nodes in column $i$ and is zero elsewhere, with box constraints $-1 \leq \alpha_i \leq 1$. EdgeGraphStars are constructed analogously for edge perturbations.

## Reachability Analysis

### Affine Operations (Exact)

GraphStar sets inherit closure under affine transformations from Star sets. All linear operations in GNN layers -- weight multiplication, adjacency normalization, bias addition, source gathering ($S$), and target scattering ($R^\top$) -- preserve the GraphStar structure **exactly**. The center and each generator matrix are transformed independently.

### ReLU Activation (Sound Over-Approximation)

The only source of over-approximation is the ReLU nonlinearity. For each neuron whose pre-activation bounds straddle zero (a *crossing neuron*), the **approx-star** method introduces linear relaxation constraints that over-approximate the ReLU, solving two LP problems per neuron for tight bounds. This avoids exponential branching while preserving soundness.

### GCN Reachability

Propagating a GraphStar through one GCN layer:

1. **Linear stage (exact)**: Apply $\tilde{D}^{-1/2} \tilde{A} \tilde{D}^{-1/2} H^{(\ell)} W^{(\ell)}$ to center and all generators
2. **ReLU stage (sound)**: Over-approximate each crossing neuron via approx-star relaxation

### GINE Reachability

GINE layers combine node and edge uncertainty in the aggregation step:

1. **Source gathering (exact)**: $S$ maps the NodeGraphStar into per-edge space
2. **Minkowski sum with edge uncertainty (exact)**: The joint expression $S[\![\Theta_X^{(\ell)}]\!] \oplus [\![\Theta_E]\!] W_e^{(\ell)}$ concatenates generators from both GraphStars, producing a set with $m + r$ generators
3. **Target scattering (exact)**: $R^\top$ sums per-edge messages at target nodes
4. **MLP with inner ReLU (sound)**: Two affine maps with approx-star ReLU relaxation
5. **Outer ReLU (sound, intermediate layers only)**: Omitted at the final layer

Since $r$ (the number of edge generators) is fixed by the number of perturbed edge features and does not grow with graph size, the additional complexity from edge uncertainty is bounded.

### Soundness Guarantee

For any GraphStar input set $\Theta_\mathcal{G}$, the computed reachable set soundly over-approximates all network outputs:

$$f_\theta(X, A, E) \in \text{Reach}(f_\theta, \Theta_\mathcal{G}) \quad \text{for all } X \in [\![\Theta_X]\!], \; E \in [\![\Theta_E]\!]$$

This holds for both GCN (node-only uncertainty) and GINE (joint node-edge uncertainty).

## Subgraph Verification

In a $K$-layer GNN with local aggregation, the output at node $v$ depends only on its **$K$-hop neighborhood**. For node-level tasks (e.g., power flow, optimal power flow), GNNV extracts the $K$-hop subgraph $\mathcal{G}_v^K$ around each target node and performs reachability on this smaller graph. This dramatically reduces the number of ReLU units encountered:

- **IEEE-118** with 3-layer GINE: a typical 3-hop subgraph contains 15--35 nodes (out of 118), reducing ReLU units by 3--8x
- Enables verification in under 90 seconds at $\epsilon = 10^{-2}$

## Safety Specifications

### Voltage Magnitude Safety (Node Regression)

For power flow tasks, predicted voltage magnitudes must remain within safe operating limits. Given perturbed inputs $X', E'$ satisfying $\|X - X'\|_\infty \leq \epsilon_{\text{node}}$ and $\|E - E'\|_\infty \leq \epsilon_{\text{edge}}$:

$$\Phi_{\text{volt},i} := \hat{V}_i(X', A, E') \in [V_{\min}, V_{\max}]$$

where $[V_{\min}, V_{\max}]$ is the permissible voltage range in per unit (p.u.), e.g., $[0.95, 1.05]$ for IEEE-24.

### Classification Robustness (Graph Classification)

For graph-level classification (e.g., cascading failure analysis), the predicted class must remain invariant under perturbation:

$$\Phi_{\text{rob}} := y'_c \geq y'_j, \quad \forall j \neq c$$

where $c$ is the nominal class and $y'$ are the perturbed logits.

## Computational Complexity

The cost is dominated by ReLU branching:

- **GCN**: $K N d_h$ ReLU units per $K$-layer network
- **GINE**: $N(K d_h + (K-1) d_{\text{out}})$ total nonlinear units (e.g., $160N$ for $K=3$, $d_h = d_{\text{out}} = 32$)

Subgraph verification reduces $N$ to the local neighborhood size, making per-node verification practical.

## Experimental Results Summary

GNNV has been evaluated on three power system tasks (PF, OPF, CFA) across IEEE-24, IEEE-39, and IEEE-118 networks, as well as graph classification benchmarks (ENZYMES, PROTEINS):

- **High verification rates**: Near-perfect robustness for OPF across all systems; 70--99% for PF depending on perturbation level
- **Edge perturbation impact minimal**: Adding edge uncertainty (1% line parameter deviation) reduces robustness by at most 0.2%, with runtime overhead of 1.1--2.9x
- **Tighter than CORA**: GNNV consistently verifies more graphs than CORA's polynomial-zonotope abstractions -- up to 21.6% more at larger perturbation budgets
- **Scalable**: Subgraph verification completes in under a second for small perturbations, and under 90 seconds for large systems at $\epsilon = 10^{-2}$

## See Also

- {doc}`/user-guide/set-representations` -- GraphStar set definition and usage
- {doc}`/user-guide/architectures` -- GCN and GINE layer reference
- {doc}`/examples/gnn` -- Worked GNN verification tutorial on IEEE bus systems
- {doc}`/application-domains` -- Power systems application domain
