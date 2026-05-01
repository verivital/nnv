# NNV 3.0 Repeatability Package — ATVA 2026 Artifact

This directory is the artifact-evaluation package for the **NNV 3.0 tool
paper** (ATVA 2026, submitted). It bundles six experiments — **FairNNV**,
**ProbVer**, **GNNV**, **VideoStar**, **ModelStar**, and **ToolComparison**
(vs the MathWorks AI Verification Library) — together with a Dockerfile
that builds a self-contained MATLAB R2024b environment and a single
`run_all.sh` driver. Each experiment folder is self-contained: bundled
models, data, vendored sources, a per-folder `README.md`, and one runner
script.

## Examples at a glance

| Folder | What it verifies | Hardware | Runtime |
|--------|------------------|----------|--------:|
| [FairNNV/](FairNNV/README.md) | Counterfactual + individual fairness on Adult-Income FC nets | CPU | ~2 min |
| [ProbVer/](ProbVer/README.md) | Probabilistic robustness of TinyYOLO via cp-star reachability | GPU + ≥48 GB RAM | ~9 min |
| [GNNV/](GNNV/README.md) | Per-node voltage-magnitude bounds on GCN/SAGE/GINE-Conv (PF/IEEE24) | CPU | ~5 min |
| [VideoStar/](VideoStar/README.md) | 3D-CNN robustness on the ZoomIn-4f video classifier | GPU | ~13 min |
| [ModelStar/](ModelStar/README.md) | Weight-perturbation reachability on an MNIST MLP | CPU | ~3 min |
| [ToolComparison/](ToolComparison/README.md) | NNV vs MathWorks AIVL on ACAS Xu p3/p4, RL, OVAL21, Collins RUL, MNIST-ResNet-8 (renders ATVA paper's Tables 5, 6, 7); AIVL optional | CPU | ~12 min (smoke) / ~3–5 h (full) |

## Quickstart

```bash
# Build (~10–12 min)
docker build -t nnv3.0 -f code/nnv/examples/NNV3.0/Dockerfile .

# Run end-to-end (~30 min on RTX 5090; needs MATLAB licence)
docker run --rm --gpus all \
    -e MLM_LICENSE_FILE=<port>@<host> \
    nnv3.0 bash run_all.sh
```

`bash run_all.sh` runs the six experiments in series, each in its own MATLAB
session (a crash in one doesn't lose the others). Per-experiment logs land
in `repeatability_logs/`; `summary.csv` records wall-clock + status.

`bash run_all.sh --help` documents env-var overrides (`NNV3_SKIP`,
`TOOLCOMPARISON_MODE=full`, `NNV3_FORCE_GPU=0/1`, …).

Copy results out:

```bash
docker cp <container>:/home/matlab/nnv/code/nnv/examples/NNV3.0 ./nnv30_results
```

## Prerequisites

- **Docker** ≥ 24 with the daemon running
- **NVIDIA driver** ≥ 535 + Container Toolkit — only for ProbVer / VideoStar
  (GNNV is CPU-only by design; without a GPU `run_all.sh` skips ProbVer
  and runs the rest)
- **MATLAB licence** the container can reach: a network server
  (`port@host` via `MLM_LICENSE_FILE`) or a `network.lic` you can mount.
  MathWorks issues 30-day eval licences that work
  (https://www.mathworks.com/campaigns/products/trials.html).
- **≥ 16 GB RAM**; **≥ 48 GB for ProbVer** (cp-star can produce ImageStars
  >31 GB; `run_all.sh` auto-skips below `NNV3_MIN_MEMORY_GB=48`,
  override with `NNV3_FORCE_MEMORY=1`)
- **~14 GB disk** for the image, ~1 GB for outputs

Host setup details (Windows WSL2 RAM cap, GPU sanity check, build-arg
`LICENSE_SERVER`) are in [§Host setup](#host-setup) at the end.

## Run individual experiments

After `docker run --gpus all -it -e MLM_LICENSE_FILE=<port>@<host> nnv3.0`
you land in `/home/matlab/nnv/code/nnv/examples/NNV3.0/`. `cd` into the
folder for the experiment you want, then run the command in its row:

| Folder | Run command (from inside the folder) | Knobs | Key outputs |
|--------|--------------------------------------|-------|-------------|
| `FairNNV` | `matlab -nodisplay -r "run('run_fairnnv.m'); exit()"` | edit script | `results/<ts>/{counterfactual,individual,timing}_*.{csv,tex}` |
| `ProbVer` | `bash run_probver.sh` | `PROBVER_NUM_SAMPLES`, `PROBVER_SEED`, `PROBVER_NRAND` | `results/results_summary.csv` |
| `GNNV` | `matlab -nodisplay -r "run_gnn_experiments(); exit()"` | `'architectures',{…}`, `'mode','node_edge'`, `'num_graphs',N`, `'node_epsilons',[…]` | `results/gnn_<ts>/{gnn_results.csv,results.mat}` |
| `VideoStar` | `matlab -nodisplay -r "run('run_zoomin_4f.m'); exit()"` | `config.sampleIndices`, `config.verAlgorithm`, `config.timeout` | `results/<ts>/eps=*.csv` |
| `ModelStar` | `matlab -nodisplay -r "run('run_expt_for_compute.m'); run('EXPT.m'); exit()"` | (paper-default sweep) | `results/MNIST_MLP.mat` + heatmap |
| `ToolComparison` | `matlab -nodisplay -r "run('run_toolcomparison.m'); exit()"` | `'mode','smoke'`, `'tools',{'nnv'}` | `acas_rl_tll/results/*.{mat,csv}`, `mnist_resnet/results/*.{mat,csv}`, `tables/out/{table_A,table_C,sanity_report}.{tex,txt}` |

ProbVer's `run_probver.sh` runs each instance in a fresh MATLAB process so
an OOM-SIGKILL on one instance doesn't kill the suite (cp-star can
overshoot 30 GB on TinyYOLO). Single-process fallback is `matlab -r
"run('run_probver.m'); exit()"`.

ToolComparison defaults to **full** when called directly (~3–5 h, renders
the ATVA paper's Tables 5, 6, and 7 plus the CAV'23 sanity report). Pass
`'mode','smoke'` for a 12 min NNV-only sanity pass on 5 ACAS + 5 MNIST
images. The AIVL Support Package is optional; without it pass
`'tools',{'nnv'}` to skip the MathWorks-side rows. Tarball recipe and
extraction details are in [`ToolComparison/README.md`](ToolComparison/README.md).

## Reference timings

Wall-clock from a Windows 11 host with RTX 5090 (32 GB, Blackwell), driver
581.95, CUDA 13. MATLAB R2024b inside the container.

| Experiment | Wall-clock | Notes |
|------------|-----------:|-------|
| FairNNV | 121 s | 100 samples × 7 ε × 2 ONNX, CPU |
| ProbVer (3 inst) | 519 s | TinyYOLO + cp-star, GPU; 3/3 unsat |
| GNNV (120 verifs) | ~5 min | 10 graphs × 3 archs × 4 ε on PF/IEEE24, CPU |
| VideoStar ZoomIn-4f | 754 s | 10 samples × 3 ε with `relax`, GPU; 7 verified / 3 unknown per ε |
| ModelStar (fc_4–fc_6) | ~3 min | MNIST MLP weight-perturbation sweep, CPU |
| ToolComparison smoke | ~12 min | 5 ACAS + 5 MNIST, NNV-only |
| ToolComparison full | ~3–5 h | NNV per-benchmark grid + AIVL, CPU |
| **End-to-end (`run_all.sh`)** | **~30 min** | All six (ToolComparison defaults to smoke inside `run_all.sh`) |

Expect 1.5–3× longer on RTX 4070 / A4000, 5–10× on CPU-only hosts.
Verdicts are hardware-independent; see
[`EXPECTED_RESULTS.md`](EXPECTED_RESULTS.md) for the cell-level
verified-counts and timing baselines.

## Host setup

### Build with licence baked in (optional)

```bash
docker build \
    -t nnv3.0 \
    -f code/nnv/examples/NNV3.0/Dockerfile \
    --build-arg LICENSE_SERVER=<port>@<host> \
    .
```

This sets `MLM_LICENSE_FILE` inside the image so you can drop the
`-e MLM_LICENSE_FILE=…` from `docker run`. The build itself does not
validate the licence; it's consumed at first MATLAB invocation.

### Windows / Docker Desktop RAM cap

WSL2 defaults to ~50% of host RAM. For ProbVer, raise it via
`%USERPROFILE%\.wslconfig`:

```ini
[wsl2]
memory=56GB
swap=32GB
processors=auto
```

Then `wsl --shutdown` in PowerShell and relaunch Docker Desktop.

### GPU sanity check

```bash
docker run --rm --gpus all nvidia/cuda:12.8.2-base-ubuntu22.04 nvidia-smi
```

If your GPU appears, `--gpus all` is wired up correctly. (Newer CUDA-13
tags also work on Blackwell hosts; pick any tag from the
[`nvidia/cuda` Docker Hub repo](https://hub.docker.com/r/nvidia/cuda/tags?name=base-ubuntu22.04).)

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `License Manager Error` on first MATLAB invocation | No licence reachable from container | `--build-arg LICENSE_SERVER=port@host` at build, or `-e MLM_LICENSE_FILE=port@host` at run, or mount `network.lic` |
| `GPU device is not supported because it has a higher compute capability…` | RTX 50-series (Blackwell, CC 12.0) is newer than MATLAB R2024b's bundled CUDA | Scripts already call `parallel.gpu.enableCUDAForwardCompatibility(true)`. On MATLAB without that API, upgrade to R2025a+ (CUDA 12) or R2025b (CUDA 13) |
| `Undefined function 'load_vnnlib'` from ProbVer | NNV path not on MATLAB path | `addpath(genpath('/home/matlab/nnv/code/nnv'))` |
| `Unrecognized function or variable 'glpk'` from ToolComparison's `lpsolver` | Host bind-mount over `code/nnv/` hides image's populated `tbxmanager/` directory | `docker cp <container>:/home/matlab/nnv/code/nnv/tbxmanager ./tbxmanager` so the bind mount has GLPK visible |
| `tbxmanager.com` fetch error during `install.m` | Mirror intermittently unreachable | Caught by the build; NNV3.0 examples don't need MPT3. Run `check_nnv_setup` to confirm core NNV |
| AIVL "Support Package required" inside ToolComparison parfeval workers | Workers don't auto-run `startup.m`; manual-extract AIVL isn't registered with the Support Package system | The verifier's MW path re-adds `aivnv` on entry. Confirm `utils/toolbox_install.m` ran and `/home/matlab/Documents/MATLAB/SupportPackages/R*/.../aivnv` exists |
| `scripts/` not found inside ToolComparison | Helpers moved to `utils/` (May 2026 refactor) | Use `cd ToolComparison/utils` |
| `npy-matlab directory not found` | Pre-built image without clone | Pull latest image, or `git clone https://github.com/kwikteam/npy-matlab code/nnv/examples/Submission/FORMALISE2025/npy-matlab` |
