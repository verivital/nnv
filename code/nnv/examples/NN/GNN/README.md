# GNN Verification (GNNV)

**Full documentation:** [GNN Tutorial](https://verivital.github.io/nnv/examples/gnn.html) · [GNN Theory](https://verivital.github.io/nnv/theory/gnn-reachability.html)

## Quick Start

```matlab
results = run_pf_verification('IEEE24', 'GCN', 0.01);
results = run_pf_verification('IEEE39', 'GINE', 0.01);
results = run_opf_verification('IEEE24', 'GINE', 0.01);
```

## Supported Architectures

- **GCN** (Graph Convolutional Network): Node-only features, uses normalized adjacency
- **GINE** (Graph Isomorphism Network with Edge features): Node + edge features

## IEEE Bus Systems

| System | Nodes | Edges |
|--------|-------|-------|
| IEEE 24-bus | 24 | 92 |
| IEEE 39-bus | 39 | 131 |
| IEEE 118-bus | 118 | 476 |

## Key Classes

`GNN`, `GCNLayer`, `GINELayer`, `GraphStar`

See the [documentation](https://verivital.github.io/nnv/examples/gnn.html) for complete walkthrough, verification status codes, and importing your own models.
