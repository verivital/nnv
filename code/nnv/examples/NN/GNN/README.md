# GNN Verification Examples

This folder contains examples demonstrating Graph Neural Network (GNN) verification using NNV. The examples use trained GNN models for power system applications: Power Flow (PF) and Optimal Power Flow (OPF) prediction.

## Directory Structure

```
GNN/
├── PowerFlow/
│   ├── run_pf_verification.m         # Helper function for all PF experiments
│   ├── verify_voltage_spec.m         # Voltage specification verification
│   ├── IEEE24/
│   │   ├── models/
│   │   │   ├── gcn_pf_ieee24.mat
│   │   │   └── gine_pf_ieee24.mat
│   │   ├── verify_pf_gcn.m
│   │   ├── verify_pf_gine.m
│   │   └── verify_pf_gine_edge_perturb.m
│   ├── IEEE39/
│   │   └── (same structure)
│   └── IEEE118/
│       └── (same structure)
├── OptimalPowerFlow/
│   ├── run_opf_verification.m        # Helper function for all OPF experiments
│   ├── IEEE24/
│   │   ├── models/
│   │   │   ├── gcn_opf_ieee24.mat
│   │   │   └── gine_opf_ieee24.mat
│   │   ├── verify_opf_gcn.m
│   │   └── verify_opf_gine.m
│   ├── IEEE39/
│   │   └── (same structure)
│   └── IEEE118/
│       └── (same structure)
└── README.md
```

## Quick Start

### Using Helper Functions (Recommended)

```matlab
% Power Flow verification
results = run_pf_verification('IEEE24', 'GCN', 0.01);
results = run_pf_verification('IEEE39', 'GINE', 0.01);

% Optimal Power Flow verification
results = run_opf_verification('IEEE24', 'GINE', 0.01);
results = run_opf_verification('IEEE118', 'GCN', 0.01);

% With edge perturbation (GINE only)
results = run_pf_verification('IEEE24', 'GINE_edge', 0.01);
```

### Running Individual Scripts

```matlab
cd('examples/NN/GNN/PowerFlow/IEEE24');
verify_pf_gcn;    % GCN verification
verify_pf_gine;   % GINE verification
```

## Model Architectures

| Architecture | Description | Layer Operation |
|--------------|-------------|-----------------|
| GCN | Graph Convolutional Network | Y = A_norm * X * W + b |
| GINE | Graph Isomorphism Network with Edge features | Y = W_node * ((1+ε)*X + Σ_j ReLU(h_j + φ_edge(e_ij))) |

Both use 3-layer architectures with ReLU activations.

## Verification Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Epsilon | 0.01 (1%) | Perturbation magnitude relative to feature range |
| Perturbed features | [1, 2] | Power injections (Pg-Pd, Qg-Qd) |
| Voltage spec | 0.95-1.05 p.u. | Safety bounds for verification |
| Reachability method | approx-star | Approximate star-based reachability |

## Verification Results

Results include verification status for each voltage node:

| Status | Code | Meaning |
|--------|------|---------|
| Verified safe | 1 | Output bounds fully within spec |
| Violated | 0 | Output bounds fully outside spec |
| Unknown (boundary) | 2 | Bounds cross specification boundary |
| Unknown (timeout) | 3 | Verification timed out |
| N/A | -1 | Non-voltage bus (not verified) |

### Example Output

```
=== GINE Power Flow Verification (IEEE24) ===
Model: .../models/gine_pf_ieee24.mat
Epsilon: 0.0100 (1.00%)
Graph: 24 nodes, 92 edges

Computing reachability...
Completed in 1.64 seconds
Samples within bounds: 10/10

=== Voltage Verification ===
Spec: 0.95 <= V <= 1.05 p.u.
  Verified safe: 7 nodes
  Violated: 0 nodes
  Unknown: 6 nodes
    - Bounds cross spec boundary: 6
    - Timeout/inconclusive: 0
```

## IEEE Bus Systems

| System | Nodes | Edges | Description |
|--------|-------|-------|-------------|
| IEEE 24 | 24 | 92 | IEEE Reliability Test System |
| IEEE 39 | 39 | 131 | New England Test System |
| IEEE 118 | 118 | 476 | IEEE 118-bus Test System |

## Importing Your Own Models

### Weight Format
NNV uses weight matrices with shape `(F_in x F_out)`:
- **PyTorch/PyG**: Uses `(F_out x F_in)` — transpose required
- **MATLAB**: Already in correct format

### Architecture Compatibility

| NNV Layer | Compatible With | Notes |
|-----------|----------------|-------|
| GCNLayer | PyTorch Geometric GCNConv | Transpose weights |
| GINELayer | Custom linear-projection GINE | NOT PyG GINEConv |

**Note**: GINELayer uses linear projections (not MLPs) for verification soundness.

## Key Classes

- `GNN`: Wrapper class for GNN networks
- `GCNLayer`: Graph Convolutional layer
- `GINELayer`: GINE layer with edge features
- `GraphStar`: Star set representation for graph node features

## References

- Kipf & Welling, "Semi-Supervised Classification with Graph Convolutional Networks", ICLR 2017
- Hu et al., "Strategies for Pre-training Graph Neural Networks", ICLR 2020
- Varbella et al., "PowerGraph: A Power Grid Benchmark Dataset for Graph Neural Networks", NIPS 2024
