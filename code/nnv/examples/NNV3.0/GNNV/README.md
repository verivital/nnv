# GNNV — Graph Neural Network Verification

Reachability-based verification of Graph Neural Networks for **AC power
flow** prediction on the IEEE 24-bus system. The script runs three
PyTorch Geometric–compatible architectures (GCN, SAGE, GINE-Conv) over
a configurable ε grid and checks per-node voltage-magnitude safety
constraints. Models load through `gnn2nnv.m`, which auto-detects the
architecture from the `.mat`'s `model_type` field and builds the right
NNV `GNN` object.

## References

- **GNNV (this work)**: Tumlin, A.M., Shao, Z., Manzanas Lopez, D.,
  Derr, T., Johnson, T.T. *Reachability-based formal verification of
  graph neural networks with node and edge features.* 2026 (under
  review).
- **Power-grid benchmark**: Varbella, A., Amara, K., Gjorgiev, B.,
  El-Assady, M., Sansavini, G. *PowerGraph: A power grid benchmark
  dataset for graph neural networks.* NeurIPS 2024 Datasets &
  Benchmarks.

## Architectures

Three GNN architectures, all with direct PyTorch Geometric analogues:

| Tag         | PyG layer  | NNV class       | Edge features used? |
|-------------|------------|-----------------|---------------------|
| `gcn`       | `GCNConv`  | `GCNLayer`      | no                  |
| `sage`      | `SAGEConv` | `SAGEConvLayer` | no (binary adjacency only) |
| `gine_conv` | `GINEConv` | `GINEConvLayer` | yes (line impedance) |

## Layout

```
NNV3.0/GNNV/
├── README.md                       This file
├── run_gnn_experiments.m           Main verification script
├── PowerFlow/IEEE24/models/
│   ├── gcn_pf_ieee24.mat
│   ├── sage_pf_ieee24.mat
│   └── gine_conv_pf_ieee24.mat
└── results/                        Timestamped output (gnn_<ts>/)
```

Only PowerFlow on IEEE 24-bus ships in this folder. The script is
extension-ready (additional grids would live at `PowerFlow/IEEE39/`,
etc.) but no other grids are included.

## Running

NNV must be installed (`code/nnv/install.m`) and on the MATLAB path
inside the Docker image.

### Default sweep (paper-aligned)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/GNNV; run_gnn_experiments"
```

Sweeps all 3 architectures × PF × IEEE24 × node-only perturbation × 4 ε
values × 10 graphs.

### Smoke (single architecture, single ε, fewer graphs)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/GNNV; \
    run_gnn_experiments('num_graphs', 3, 'architectures', {'gcn'}, \
                        'node_epsilons', [1e-3])"
```

## Configuration parameters

| Parameter        | Default                       | Notes |
|------------------|-------------------------------|-------|
| `architectures`  | `{'gcn','sage','gine_conv'}`  | Subset to run |
| `grid`           | `'ieee24'`                    | Only IEEE24 ships |
| `task`           | `'pf'`                        | Only PowerFlow ships |
| `mode`           | `'node_only'`                 | `'node_only'` or `'node_edge'` (GINE-Conv + PF only) |
| `num_graphs`     | `10`                          | Pre-filtered to "safe" graphs (GT voltages within voltage spec) |
| `node_epsilons`  | `[1e-5, 1e-4, 1e-3, 1e-2]`    | Relative perturbation on (P, Q) |
| `edge_epsilons`  | `[1e-3, 1e-2]`                | Used only when `mode='node_edge'` |
| `parallel`       | `false`                       | Set `true` for `parfor` |
| `num_workers`    | `0` (auto)                    | When `parallel=true` |

## Outputs

A timestamped subfolder `results/gnn_<yymmdd-HHMMSS>/` is created per
run. Inside:

- `gnn_results.csv` — flat row-per-config CSV. Columns: Architecture,
  Task, Mode, Grid, Node_Epsilon, Edge_Epsilon, Total_Graphs,
  Safe_Graphs, Skipped_Graphs, Avg_Time_s, Total_Verified,
  Total_Unknown, Total_Violated, Total_Voltage_Nodes, Pct_Verified,
  Mean_Bound_Width, Max_Bound_Width
- `results.mat` — full nested struct keyed `arch.task.mode.grid.eps_*`
- `experiments.log` — MATLAB diary

Plots and LaTeX tables for the paper are produced post-hoc from `gnn_results.csv` with external tools.

## Expected runtime

Measured on an NVIDIA RTX 4060 host running the `nnv3.0` Docker image
(MATLAB R2024b), CPU verification:

- **Smoke** (1 arch × 1 ε × 3 graphs): ~10 s
- **Default sweep** (3 archs × 4 ε × 10 graphs, PF only): **~3.5 minutes**
