# FairNNV — Fairness Verification of the Adult-Income Classifier

Exact reachability-based verification of two notions of fairness on
binary classifiers trained on the UCI Adult-Income dataset:

- **Counterfactual fairness** — flipping a sensitive attribute (e.g.
  sex) must not change the prediction. Verified at ε = 0.
- **Individual fairness** — for every input, no perturbation within an
  ε-ball (combined with a flip of the sensitive attribute) changes the
  prediction. Verified across multiple ε values.

Each verdict is summarized as a Verified Fairness (VF) score: the
proportion of test samples for which fairness is formally certified.

## References

- **FairNNV (this work)**: Tumlin, A.M., Manzanas Lopez, D., Robinette,
  P., Zhao, Y., Derr, T., Johnson, T.T. *FairNNV: The neural network
  verification tool for certifying fairness.* Proceedings of the 5th
  ACM International Conference on AI in Finance (ICAIF '24), 2024.
- **Counterfactual fairness definition**: Kusner, M.J., Loftus, J.R.,
  Russell, C., Silva, R. *Counterfactual fairness.* NeurIPS 2017.
- **Adult-Income dataset**: Dheeru & Efi. *UCI Machine Learning
  Repository — Adult.* 2017.

## Models

Two ONNX classifiers in `models/`:

| Model | Architecture           | Notes                       |
|-------|------------------------|-----------------------------|
| AC-1  | 13 → 16 → 8 → 2        | "Small": two narrow hidden layers |
| AC-3  | 13 → 50 → 2            | "Medium": one wider hidden layer  |

## Layout

```
NNV3.0/FairNNV/
├── README.md
├── run_fairnnv.m       Top-level runner; sets config and chains the steps
├── adult_verify.m      Loads ONNX, runs reachability + verification, writes CSVs
├── plot_results.m      Reads the latest CSVs, generates figures + LaTeX tables
├── models/
│   ├── AC-1.onnx
│   └── AC-3.onnx
├── data/
│   └── adult_data.mat
└── results/            Timestamped output (<yymmdd-HHMMSS>/)
```

`adult_verify.m` and `plot_results.m` can also run standalone — they
fall back to default paths in this folder when `config` is not already
in scope.

## Running

NNV must be installed (`code/nnv/install.m`) and on the MATLAB path
inside the Docker image. The runner adds NNV automatically using a
path relative to this folder.

### Default sweep

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/FairNNV; run_fairnnv"
```

Verifies AC-1 and AC-3 on 100 observations, counterfactual fairness
(ε = 0) plus individual fairness across the paper's ε grid.

### Smoke

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/FairNNV; \
    config.modelList = {'AC-1'}; config.numObs = 10; \
    config.epsilon_individual = [0.01]; config.timeout = 120; \
    run_fairnnv"
```

## Configuration parameters

Edit the `CONFIGURATION` block at the top of [`run_fairnnv.m`](run_fairnnv.m),
or pre-set fields in the `config` struct before calling the runner
(the runner's config-guard preserves any caller-supplied values):

| Field                          | Default                                | Effect |
|--------------------------------|----------------------------------------|--------|
| `config.modelList`             | `{'AC-1', 'AC-3'}`                     | Which models to verify |
| `config.numObs`                | `100`                                  | Number of test observations |
| `config.randomSeed`            | `500`                                  | RNG seed |
| `config.timeout`               | `600`                                  | Per-epsilon timeout (s) |
| `config.epsilon_counterfactual`| `[0.0]`                                | ε grid for counterfactual |
| `config.epsilon_individual`    | `[0.01,0.02,0.03,0.05,0.07,0.1]`       | ε grid for individual |
| `config.savePNG/savePDF`       | `true`                                 | Figure formats |

## Outputs

A timestamped subfolder `results/<yymmdd-HHMMSS>/` is created per run
and contains:

- `counterfactual_<ts>.csv` — per-model fair / unfair %
- `individual_<ts>.csv`     — per-model × ε fair / unfair / unknown %
- `timing_<ts>.csv`         — per-model × ε total + per-sample time
- `counterfactual_table.tex` — booktabs-style LaTeX table
- `individual_fairness_combined.png` / `.pdf` — area plot across models
- `timing_table.tex`         — LaTeX timing table

## Expected runtime

Measured on an NVIDIA RTX 4060 host running the `nnv3.0` Docker image
(MATLAB R2024b), CPU verification:

- **Smoke** (AC-1 only, 10 obs, ε ∈ {0.01}): **~30 s**
- **Default sweep** (2 models, 7 ε values, 100 obs, including plotting): **~2 minutes**
