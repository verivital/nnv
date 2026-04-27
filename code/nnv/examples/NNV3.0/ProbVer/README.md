# ProbVer — Probabilistic Verification of TinyYOLO

Probabilistic, model-agnostic, data-driven verification of the
**TinyYOLO** object detector against the YOLO-2023 VNN-COMP benchmark.
The pipeline combines randomized falsification with conformal-inference
reachability via NNV's **cp-star** method. For each VNNLIB property:

1. Random-sample inputs to find a counterexample (SAT short-circuit).
2. If no counterexample, train a surrogate linear model from samples,
   compute conformal bounds on the output, and check whether the
   bounds satisfy the specification (UNSAT → property holds with the
   declared coverage / confidence).

The approach scales to networks where exact reachability is
intractable, providing probabilistic guarantees instead of soundness.

## References

- **Probabilistic verification (this work)**: Hashemi, N., Sasaki, S.,
  Oguz, I., Ma, M., Johnson, T.T. *Scaling data-driven probabilistic
  robustness analysis for semantic segmentation neural networks.*
  NeurIPS 2025.
- **Benchmark**: Brix, C., Bak, S., Liu, C., Johnson, T.T. *The fourth
  international verification of neural networks competition (VNN-COMP
  2023): summary and results.* arXiv:2312.16760, 2023.

## Layout

```
NNV3.0/ProbVer/
├── README.md
├── run_probver.m                  Main verification script
└── yolo_2023/                     Self-contained benchmark
    ├── instances.csv              72 instances (onnx, vnnlib, timeout)
    ├── onnx/TinyYOLO.onnx.gz      Compressed model
    └── vnnlib/*.vnnlib.gz         72 compressed properties
```

`run_probver.m` decompresses the `.gz` files on first run.

## Running

NNV must be installed (`code/nnv/install.m`) and on the MATLAB path
inside the Docker image. The Python helpers used by cp-star
(`Direction_trainer.py`, `Trainer_Linear.py`, `Trainer_ReLU.py`)
require the host `.venv` to exist; see *One-time host setup* in the
top-level [`NNV3.0/README.md`](../README.md).

### Default sweep (paper-aligned)

```bash
docker run --rm --gpus all \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv/code/nnv/examples/NNV3.0/ProbVer \
    nnv3.0 matlab -nodisplay -batch "run_probver"
```

3 instances are randomly selected (seed 42) from the 72-instance
benchmark; each is verified at ε = 1/255.

### Smoke (single instance, fewer falsification samples)

```bash
docker run --rm --gpus all \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv/code/nnv/examples/NNV3.0/ProbVer \
    nnv3.0 matlab -nodisplay -batch \
    "numSamples = 1; nRand = 50; run_probver"
```

## Configuration parameters

Edit the `CONFIGURATION` block at the top of
[`run_probver.m`](run_probver.m), or set the variables before invoking
the runner (the runner has a config-guard for these):

| Variable      | Default | Notes |
|---------------|---------|-------|
| `numSamples`  | `3`     | Number of instances to verify (out of 72) |
| `randomSeed`  | `42`    | RNG seed for instance sampling            |
| `nRand`       | `100`   | Random samples per instance for falsification |

cp-star reachability options (set inside `run_probver.m` `reachOptions`):

| Option              | Default      | Meaning |
|---------------------|--------------|---------|
| `train_epochs`      | `200`        | Surrogate linear model training epochs |
| `train_lr`          | `0.0001`     | Learning rate |
| `coverage`          | `0.999`      | Probabilistic coverage |
| `confidence`        | `0.999`      | Probabilistic confidence |
| `train_mode`        | `'Linear'`   | Surrogate type (`'Linear'` or `'ReLU'`) |
| `surrogate_dim`     | `[-1, -1]`   | Auto-pick from input shape |
| `threshold_normal`  | `1e-5`       | Normality test threshold |
| `reachMethod`       | `'cp-star'`  | NNV reach method |

## Outputs

A timestamped subfolder `results/<yymmdd-HHMMSS>/` is created per run.
Inside:

- `results_summary.csv` — one row per verified instance:
  `Index, ONNX, VNNLIB, Status, StatusCode, Time(s), Error`. Status
  codes: `0` = SAT (counterexample), `1` = UNSAT (verified), `2` =
  unknown, `-1` = error.

## Expected runtime

Measured on an NVIDIA RTX 4060 host running the `nnv3.0` Docker image
(MATLAB R2024b). The torch wheels installed in the host `.venv`
require CUDA 13; with this driver the GPU detection fails and torch
falls back to CPU — still runs, just slower.

- **Smoke** (1 instance, `nRand=50`, CPU surrogate training): **~3 minutes**
- **Default sweep** (3 instances, `nRand=100`, CPU): **~8 minutes**

