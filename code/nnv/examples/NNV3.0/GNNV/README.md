# GNN Verification (NNV3.0/GNNV)

Reachability-based verification of Graph Neural Networks for power-flow
prediction on IEEE bus systems. Models are PyTorch Geometric–compatible
architectures, exported to NNV via the `gnn2nnv` loader.

## Architectures

Three GNN architectures with direct PyG analogues:

| Tag         | PyG layer    | NNV class         | Edge features? |
|-------------|--------------|-------------------|----------------|
| `gcn`       | `GCNConv`    | `GCNLayer`        | no             |
| `sage`      | `SAGEConv`   | `SAGEConvLayer`   | no             |
| `gine_conv` | `GINEConv`   | `GINEConvLayer`   | yes            |

Models trained externally and shipped here as `.mat` files. The
`model_type` field inside each `.mat` is what `gnn2nnv` uses to dispatch
to the right NNV layer construction.

## Tasks and grids

- **Tasks**: PowerFlow (`pf`), OptimalPowerFlow (`opf`).
- **Grid (shipped)**: IEEE 24-bus only. Folder layout
  (`PowerFlow/IEEE39/`, `IEEE118/`) is reserved for future grids — drop
  more `.mat` files in and the harness picks them up via
  `'grid','ieee39'`.

## Layout

```
NNV3.0/GNNV/
├── run_gnn_experiments.m        Unified verification harness
├── soundness_check.m            MC-sampling vs reachable bounds
├── PowerFlow/IEEE24/models/     gcn|sage|gine_conv _pf_ieee24.mat
├── OptimalPowerFlow/IEEE24/models/  gcn|sage|gine_conv _opf_ieee24.mat
├── results/                     Timestamped output directories
└── figures/                     Reserved for post-hoc plots
```

## Running

NNV must be installed (`code/nnv/install.m`) and on the MATLAB path.

### Default sweep (all archs × both tasks × IEEE24)

```matlab
run_gnn_experiments();
```

Outputs go to `results/gnn_<timestamp>/`:
- `results.mat` — full nested struct keyed `arch.task.mode.grid.eps_*`
- `gnn_results.csv` — flat row-per-config CSV
- `experiments.log` — diary

### Smaller smoke test

```matlab
run_gnn_experiments('num_graphs', 5, ...
                    'architectures', {'gcn'}, ...
                    'task', 'pf', ...
                    'mode', 'node_only', ...
                    'node_epsilons', [1e-3]);
```

### Edge perturbation (PF only, GINE-Conv only)

```matlab
run_gnn_experiments('architectures', {'gine_conv'}, ...
                    'task', 'pf', ...
                    'mode', 'node_edge', ...
                    'edge_epsilons', [1e-3, 5e-3]);
```

### Configuration parameters

| Parameter        | Default                       | Notes |
|------------------|-------------------------------|-------|
| `architectures`  | `{'gcn','sage','gine_conv'}`  | Subset to run |
| `grid`           | `'ieee24'`                    | Or `'all'` if more grids ship later |
| `task`           | `'all'`                       | Or `'pf'` / `'opf'` |
| `mode`           | `'all'`                       | `'node_only'` or `'node_edge'` |
| `num_graphs`     | `100`                         | Pre-filtered to "safe" graphs (GT voltages within spec) |
| `node_epsilons`  | `[1e-5,1e-4,1e-3,1e-2]`       | Relative perturbation on (P, Q) |
| `edge_epsilons`  | `[1e-3,1e-2]`                 | Relative perturbation on impedance, GINE-Conv only |
| `parallel`       | `false`                       | Set `true` for `parfor` |
| `num_workers`    | `0` (auto)                    | Workers when `parallel=true` |

## Soundness check

Validates that reachable bounds over-approximate the true output:
```matlab
soundness_check();                   % default: 3 graphs × 50 MC samples
soundness_check('archs', {'gcn'});   % single arch
```

For each shipped `(arch, task, grid)` config, the script:
1. Verifies nominal accuracy: `gnn.evaluate(X)` matches stored `Y_test`.
2. Runs reachability for ε = 0.005 and Monte Carlo samples 50 random
   perturbations within the ε-ball — every sample's voltage output
   must lie inside the reachable bounds.
3. Reports mean / max bound widths.

## Plotting / tables

`gnn_results.csv` is the single source of truth for results. Plots and
LaTeX tables are intentionally not generated here — produce them with
your own tooling (Python / pandas / matplotlib / siunitx) from the CSV.

## Adding new architectures or grids

1. Train a model in PyTorch Geometric. Export to `.mat` via the
   `mat_exporter.py` / `weight_converter.py` pipeline (lives in the
   external research repo, not in NNV).
2. Drop the `.mat` into the right folder
   (`<Task>/<GRID>/models/<arch>_<task>_<grid>.mat`) and ensure
   `model_type` is set inside the `.mat`.
3. If the architecture is not already supported by `gnn2nnv`, extend
   `code/nnv/engine/utils/gnn2nnv.m` and add a corresponding
   `<Arch>Layer.m` under `code/nnv/engine/nn/layers/`.
4. Add the architecture tag to the `architectures` arg of
   `run_gnn_experiments`.
