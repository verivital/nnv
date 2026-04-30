# NNV 3.0 Repeatability Package

This directory bundles the five experiments demonstrating NNV 3.0's new
capabilities — **FairNNV**, **ProbVer**, **GNNV**, **VideoStar**, and
**ModelStar** — together with a Dockerfile that builds a self-contained
MATLAB R2024b environment and a single `run_all.sh` driver that executes
all five end-to-end. Each example folder is self-contained: bundled
models, data, vendored sources, a per-folder `README.md`, and a single
top-level runner script.

## Examples at a glance

| Folder      | What it verifies                                              | Hardware     |
|-------------|---------------------------------------------------------------|--------------|
| [FairNNV/](FairNNV/README.md)       | Counterfactual + individual fairness on Adult-Income FC nets         | CPU          |
| [ProbVer/](ProbVer/README.md)       | Probabilistic robustness of TinyYOLO via cp-star reachability        | GPU required |
| [GNNV/](GNNV/README.md)             | Per-node voltage-magnitude bounds on GCN/SAGE/GINE-Conv (PF/IEEE24)  | CPU/GPU      |
| [VideoStar/](VideoStar/README.md)   | 3D-CNN robustness on the ZoomIn-4f video classifier                  | GPU required |
| [ModelStar/](ModelStar/README.md)   | Weight-perturbation reachability on an MNIST MLP                     | CPU          |
| [ToolComparison/](ToolComparison/README.md) | NNV vs MathWorks AI Verification Toolbox on ACAS Xu / RL / TLL / MNIST-ResNet-8 | CPU; AI Verification Toolbox required |

## Prerequisites

Before the first `docker build`, the host must already have:

- **Docker** ≥ 24, with the daemon running
- **NVIDIA driver** ≥ 535 (CUDA 12+) on the host (only required for the GPU
  experiments — ProbVer, GNNV, VideoStar)
- **NVIDIA Container Toolkit** so `docker run --gpus all` works. On Windows,
  this means **WSL2** with at least one registered Linux distro (e.g. Ubuntu)
  and Docker Desktop's WSL Integration enabled
- **MATLAB licence** that the container can reach. Either a network licence
  server (`port@host`) or a node-locked licence file you can mount into the
  container. The Dockerfile takes a `LICENSE_SERVER` build arg; if you set
  it, it's baked in as `MLM_LICENSE_FILE`. Otherwise, supply
  `-e MLM_LICENSE_FILE=...` (or mount `network.lic`) at `docker run` time.
- **Disk**: ~14 GB for the image, plus ~1 GB for experiment artefacts.
- **RAM**: ≥ 16 GB to run FairNNV / GNNV / VideoStar / ModelStar; **≥ 48 GB recommended for ProbVer**.
  ProbVer's cp-star reachability on TinyYOLO can produce intermediate ImageStars
  whose peak size exceeds 31 GB. With less than that, MATLAB is liable to be
  SIGKILLed mid-run on certain instances. `run_all.sh` measures container RAM at
  startup and auto-skips ProbVer below the threshold (override with
  `NNV3_FORCE_MEMORY=1`). Tune the threshold with `NNV3_MIN_MEMORY_GB`.

  **Windows / Docker Desktop users**: by default WSL2 gives the container only
  ~50% of host RAM. Raise the cap by creating `%USERPROFILE%\.wslconfig`:

  ```ini
  [wsl2]
  memory=56GB
  swap=32GB
  processors=auto
  ```

  Then run `wsl --shutdown` in PowerShell and relaunch Docker Desktop.

Quick GPU sanity check:

```bash
docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi
```

If your GPU appears, `--gpus all` is wired up correctly.

## Build

From the **repository root** (so `.dockerignore` and the entire NNV checkout
are available to the build):

```bash
docker build \
    -t nnv3.0 \
    -f code/nnv/examples/NNV3.0/Dockerfile \
    --build-arg LICENSE_SERVER=<port>@<host> \
    .
```

Omit `--build-arg LICENSE_SERVER=...` if you'd rather provide the licence at
run time.

The build is roughly:

| Step                              | Wall-clock (RTX 5090, 1 Gbps link) |
|-----------------------------------|------------------------------------|
| Context copy + base image fetch   | ~1 min  (with `.dockerignore` — 5 min without) |
| `mpm install` MATLAB toolboxes    | ~3 min  |
| pip install Torch + CUDA wheels   | ~3 min  |
| `matlab -batch install.m`         | ~1 min  |
| Image export                      | ~3 min  |
| **Total**                         | **~10–12 min** |

## Run

```bash
docker run --gpus all -it nnv3.0
```

You'll land in `/home/matlab/nnv/code/nnv/examples/NNV3.0` with NNV already on
the MATLAB path.

### One-shot: everything

```bash
bash run_all.sh
```

Runs FairNNV, then ProbVer, then GNNV, then VideoStar, then ModelStar, each
in its own MATLAB session. Per-experiment logs land in `repeatability_logs/`,
and a final `summary.csv` records wall-clock time and exit status per
experiment.

Skip individual experiments with `NNV3_SKIP="probver videostar" bash run_all.sh`.

**GPU autodetection.** `run_all.sh` checks for an NVIDIA GPU via `nvidia-smi -L`
at startup and adapts:

- **GPU present** — runs the full suite (the README's reference timings).
- **No GPU** — prints a "deferring to CPU" notice and:
  - Skips **ProbVer** (cp-star reachability trains a surrogate network on
    CUDA via the Python venv; there is no CPU fallback for it). It appears
    in `summary.csv` as `skipped,0,no GPU (CPU fallback unsupported for
    cp-star)`.
  - Runs **FairNNV** and **ModelStar** (CPU-only by design — no change).
  - Runs **GNNV** and **VideoStar** as CPU best-effort. NNV's reachability
    is mostly CPU-side, but expect substantially longer wall-clock than
    the GPU baselines below.

Override the autodetection with `NNV3_FORCE_GPU=0` (force-skip ProbVer even
on a GPU host) or `NNV3_FORCE_GPU=1` (assume a GPU exists and try ProbVer
anyway).

**Memory autodetection.** `run_all.sh` reads `/proc/meminfo` at startup and
auto-skips ProbVer below `NNV3_MIN_MEMORY_GB` (default **48 GB**) — see the
RAM bullet under **Prerequisites**. Override with `NNV3_FORCE_MEMORY=1` to
attempt ProbVer regardless (you'll see `oom` rows in `results_summary.csv`
for instances that exceed the available memory).

### Individual experiments

#### FairNNV

```bash
cd FairNNV
matlab -nodisplay -r "run('run_fairnnv.m'); exit()"
```

Verifies counterfactual + individual fairness on the Adult-Income tabular
classifier across two ONNX models (`AC-1`, `AC-3`), 100 samples, and 7
ε values. Bundled models live in `FairNNV/models/`, the dataset in
`FairNNV/data/adult_data.mat`. CPU-only.

Outputs (in `FairNNV/results/<timestamp>/`):
- `counterfactual_*.csv` — counterfactual fairness results
- `individual_*.csv` — individual fairness results across ε
- `timing_*.csv`, `timing_table.tex`, `counterfactual_table.tex`
- `individual_fairness_combined.{png,pdf}`

#### ProbVer (requires `--gpus all`)

```bash
cd ProbVer
bash run_probver.sh
```

The bash driver verifies a random subset (default `PROBVER_NUM_SAMPLES = 3`)
of the TinyYOLO `yolo_2023` benchmark using the cp-star reachability method.
**Each instance runs in its own MATLAB process**, so if `Prob_reach` overflows
host memory and the kernel SIGKILLs MATLAB on one instance (e.g. the
`TinyYOLO_prop_000277_eps_1_255` ImageStar that warns `"large for memory"`),
the next instance starts in a fresh address space and the run continues. The
crashed instance is recorded in `results_summary.csv` with status `oom`.

Tunables (env vars):
- `PROBVER_NUM_SAMPLES` — number of instances (default 3)
- `PROBVER_SEED` — instance-selection seed (default 42)
- `PROBVER_NRAND` — falsification samples (default 100)

The single-MATLAB legacy path (`matlab -nodisplay -r "run('run_probver.m');
exit()"`) still works and is useful for debugging a single instance, but it
shares MATLAB state across instances and is more vulnerable to mid-run OOMs.

#### GNNV

```bash
cd GNNV
matlab -nodisplay -r "run_gnn_experiments(); exit()"
```

Verifies per-node voltage-magnitude bounds on the AC power-flow task for
the IEEE 24-bus grid, across three PyTorch Geometric–compatible
architectures (`GCNConv`, `SAGEConv`, `GINEConv`) and 4 ε perturbation
levels on 10 test graphs (120 verifications by default). Models load via
`gnn2nnv.m`, which auto-detects the architecture from the `.mat`'s
`model_type` field. Bundled models live in
`GNNV/PowerFlow/IEEE24/models/`.

Outputs land in `GNNV/results/gnn_<timestamp>/`:
- `gnn_results.csv` — flat per-row CSV (arch × task × eps × graph)
- `results.mat` — full nested results struct
- `experiments.log` — diary

Options (positional kwargs to `run_gnn_experiments`):
- `'architectures', {'gcn','sage'}` — subset of architectures
- `'mode', 'node_edge'` — opt-in edge-feature perturbation (GINE-Conv only)
- `'num_graphs', 5` — smoke
- `'node_epsilons', [1e-3, 5e-3]` — custom ε grid

#### VideoStar ZoomIn-4f (requires `--gpus all`)

```bash
cd VideoStar
matlab -nodisplay -r "run('run_zoomin_4f.m'); exit()"
```

Verifies the bundled `zoomin_4f.onnx` 3D-CNN video classifier on the first
10 ZoomIn test samples across ε ∈ {1/255, 2/255, 3/255}. All assets sit
next to the script: model in `VideoStar/models/`, frames in
`VideoStar/data/ZoomIn/`, vendored vvn helpers in `VideoStar/src/vvn/`,
and `VideoStar/npy-matlab/` for `.npy` reads. Configuration sits in the
script's top-level `config` struct: change `config.sampleIndices`,
`config.verAlgorithm` (`'relax'` or `'approx'`), or `config.timeout` as
needed.

Results: `VideoStar/results/<timestamp>/eps=*.csv`.

#### ModelStar

```bash
cd ModelStar
matlab -nodisplay -r "run('run_expt_for_compute.m'); run('EXPT.m'); exit()"
```

Reproduces the weight-perturbation reachability experiment on the MNIST
MLP: each fully-connected layer is perturbed independently across a sweep
of magnitudes, and per-layer robustness is verified using a built-in
weight-perturbation reachability operator. The runner reads the bundled
`mnist_model_fc.mat`, builds an experiment template programmatically (no
external YAML package required), and writes per-layer .mat results;
`EXPT.m` then renders the paper's heatmap.

Outputs:
- `ModelStar/results/MNIST_MLP.mat` — per-fraction verification results
- `ModelStar/runtime/` — per-experiment runtime CSVs
- Heatmap printed by `EXPT.m`

#### ToolComparison

```bash
cd ToolComparison
matlab -nodisplay -r "run('run_toolcomparison.m'); exit()"
```

Head-to-head between NNV's star-set reachability and the MathWorks AI
Verification Toolbox on two regimes: feed-forward networks with VNNLIB
half-space specs (ACAS Xu / RL / TLLverify) and a residual classifier
(MNIST-ResNet-8). Algorithm grid mirrors NNV 2.0 / CAV'23 — `approx-star`,
`relax-star-{range,area}-{25,50,75,100}`, `exact-star`. Renders Table A
(FC-net half) and Table C (ResNet half) plus a CAV'23 cross-check.

The MW-side comparison requires the **AI Verification Toolbox** Support
Package; see `ToolComparison/README.md` for the tarball recipe and how
to stage it for `scripts/toolbox_install.m` to extract. NNV-only runs
work without it.

```matlab
% Smoke (~10 min, NNV-only):
matlab -nodisplay -r "cd ToolComparison; run_toolcomparison('mode','smoke'); exit()"
```

Outputs:
- `ToolComparison/acas_rl_tll/results/results_<benchmark>.{mat,csv}`
- `ToolComparison/mnist_resnet/results/expC_<model>.{mat,csv}`
- `ToolComparison/tables/out/{table_A,table_C}.{tex,txt}`,
  `ToolComparison/tables/out/sanity_report.txt`

## Copy results out of the container

```bash
docker cp <container_id>:/home/matlab/nnv/code/nnv/examples/NNV3.0 ./nnv30_results
```

(`docker ps` to find the container ID.)

## Reference timings on this machine

The table below records wall-clock times measured on a Windows 11 host with an
RTX 5090 (32 GB VRAM, Blackwell, CC 12.0), driver 581.95, CUDA 13. MATLAB R2024b
runs inside the container; the same host has MATLAB R2025b natively but it
isn't used for these numbers.

| Experiment              | Wall-clock | Notes                                                                                                |
|-------------------------|-----------:|------------------------------------------------------------------------------------------------------|
| FairNNV                 |     121 s  | 100 samples × 7 ε × 2 ONNX models. CPU only.                                                         |
| ProbVer (3 instances)   |     519 s  | TinyYOLO + cp-star reach (per-instance MATLAB process). All 3 instances `unsat` (145.1, 148.1, 144.6 s). GPU. |
| GNNV (120 verifs)       |     ~5 min | 10 graphs × 3 architectures × 4 ε on PF/IEEE24. Mostly CPU-bound under MATLAB's reach.               |
| VideoStar ZoomIn-4f     |     754 s  | 10 samples × 3 ε with `relax` algorithm. GPU. 7 verified / 0 violated / 3 unknown per ε.             |
| ModelStar (fc_4–fc_6)   |     ~3 min | Last-three-layer weight-perturbation sweep on MNIST MLP. CPU only.                                   |
| ToolComparison smoke    |     ~12 min| 5 ACAS networks + 5 MNIST images, NNV-only (`approx-star` + `relax-star-*-50`). CPU only.            |
| ToolComparison full     |     ~6 h   | NNV grid (`approx-star` + 4 relax-star factors + `exact-star`) + AIVL on both halves. CPU only.      |
| **End-to-end suite**    |  **~30 min**| All five via `run_all.sh`, RTX 5090 host (excludes ToolComparison full).                            |

These numbers are intended as a baseline. Expect roughly 1.5–3× longer on a
mid-range workstation GPU (RTX 4070 / A4000) and 5–10× longer on CPU-only hosts.

## Troubleshooting

**`License Manager Error` on first `matlab` invocation.** The build doesn't
validate the licence; it's consumed at first MATLAB run. Pass
`--build-arg LICENSE_SERVER=port@host` at build time, or
`-e MLM_LICENSE_FILE=port@host` at run time, or mount a `network.lic`.

**`GPU device is not supported because it has a higher compute capability...`**
RTX 50-series (Blackwell, CC 12.0) is too new for the CUDA libraries shipped
with MATLAB R2024b. The experiment scripts call
`parallel.gpu.enableCUDAForwardCompatibility(true)` automatically. If you
are running a different MATLAB version that lacks this API, upgrade to
MATLAB R2025a+ (CUDA 12+) or R2025b (CUDA 13).

**`Undefined function 'load_vnnlib'` from ProbVer.** The script auto-bootstraps
the NNV path from its own location, so this should not happen if launched
from the `ProbVer/` directory. If it does, run
`addpath(genpath('/home/matlab/nnv/code/nnv'))` first.

**`tbxmanager.com` URL fetch error during `install.m`.** That mirror is
sometimes unreachable. The build catches it; modern NNV3.0 examples do not
require MPT3, so the install completes anyway. Run `check_nnv_setup` to
confirm core NNV is healthy before invoking experiments.

**`npy-matlab directory not found`.** Pull the latest image — the Dockerfile
now clones `kwikteam/npy-matlab` during the build. For pre-built images,
manually clone it into
`code/nnv/examples/Submission/FORMALISE2025/npy-matlab/`.
