# GNN Verification Examples

This folder contains examples demonstrating Graph Neural Network (GNN) verification using NNV. The examples use trained GNN models for power system applications: Power Flow (PF) and Optimal Power Flow (OPF) prediction.

## Directory Structure

```
GNN/
├── PowerFlow/
│   ├── verify_voltage_spec.m      # Shared voltage specification helper
│   ├── IEEE24/
│   │   ├── models/
│   │   │   ├── gcn_pf_ieee24_run4_seed131.mat
│   │   │   └── gine_edgelist_pf_ieee24_run4_seed131.mat
│   │   ├── verify_pf_gcn.m
│   │   └── verify_pf_gine.m
│   ├── IEEE39/
│   │   ├── models/
│   │   ├── verify_pf_gcn.m
│   │   └── verify_pf_gine.m
│   └── IEEE118/
│       ├── models/
│       ├── verify_pf_gcn.m
│       └── verify_pf_gine.m
├── OptimalPowerFlow/
│   ├── IEEE24/
│   ├── IEEE39/
│   └── IEEE118/
└── README.md
```

## Model Architectures

### GCN (Graph Convolutional Network)
- 3-layer GCN architecture
- Uses normalized adjacency matrix for message passing
- Layer operation: Y = A_norm * X * W + b

### GINE (Graph Isomorphism Network with Edge features)
- 3-layer GINE architecture with edge-list representation
- Incorporates edge features in message passing
- Layer operation: Y = W_node * ((1+ε)*X + Σ_j ReLU(h_j + φ_edge(e_ij)))


## IEEE Bus Systems

The power system datasets used in this work are drawn from the **PowerGraph** benchmark suite [Varbella et al., 2024], which provides standardized graph-based representations of electrical power networks for graph machine learning and verification tasks.

| System    | Nodes | Description                         |
|-----------|-------|-------------------------------------|
| IEEE 24   | 24    | IEEE Reliability Test System        |
| IEEE 39   | 39    | New England Test System             |
| IEEE 118  | 118   | IEEE 118-bus Test System            |

## Running Examples

1. Start MATLAB and navigate to the NNV root directory:
```matlab
cd('/path/to/nnv/code/nnv');
startup_nnv;
```

2. Navigate to an example folder and run:
```matlab
cd('examples/NN/GNN/PowerFlow/IEEE24');
verify_pf_gcn;   % GCN verification
verify_pf_gine;  % GINE verification
```

## Verification Workflow

Each example script performs:
1. Loads the trained model and extracts weights
2. Creates GCN or GINE layers
3. Creates a GNN wrapper with graph structure
4. Defines input perturbation (1% of feature range)
5. Creates a GraphStar input set
6. Computes reachability using approx-star method
7. Reports output bound statistics
8. Validates that random samples are within computed bounds
9. **Verifies voltage magnitude specification (0.95-1.05 p.u.)**

## Voltage Specification Verification

The examples verify that predicted voltage magnitudes satisfy grid safety constraints:

- **Specification**: 0.95 ≤ V ≤ 1.05 per-unit
- **Method**: Output bounds are checked against normalized specification bounds
- **Results**:
  - `Verified safe`: All possible outputs within spec
  - `Violated`: All possible outputs outside spec
  - `Unknown`: Output bounds cross spec boundary
  - `N/A`: Node does not predict voltage (non-voltage bus)

## Example Output

```
=== GCN Power Flow Verification (IEEE 24-bus) ===
Layer dimensions:
  Layer 1: 8 -> 64
  Layer 2: 64 -> 64
  Layer 3: 64 -> 4

Graph structure: 24 nodes

=== Computing Reachability ===
Completed in 0.1234 seconds

Output bounds - Mean: 0.012345, Max: 0.045678
Samples within bounds: 10/10

=== Voltage Specification Verification ===
Specification: 0.95 <= V <= 1.05 p.u.
  Verified safe: 18 nodes
  Violated: 0 nodes
  Unknown: 2 nodes
  N/A (non-voltage bus): 4 nodes

=== Complete ===
```

## Key Classes Used

- `GNN`: Wrapper class for GNN networks
- `GCNLayer`: Graph Convolutional layer
- `GINELayer`: GINE layer with edge features
- `GraphStar`: Star set representation for graph node features

## Helper Functions

- `verify_voltage_spec.m`: Verifies voltage magnitude bounds on GNN output
  - Inputs: GraphStar output, model data, v_min, v_max
  - Returns: Per-node verification result array

## References

- Kipf & Welling, "Semi-Supervised Classification with Graph Convolutional Networks", ICLR 2017
- Hu et al., "Strategies for Pre-training Graph Neural Networks", ICLR 2020
- Varbella et al., "PowerGraph: A Power Grid Benchmark Dataset for Graph Neural Networks", NIPS 2024