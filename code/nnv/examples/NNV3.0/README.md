# NNV3 — Tool Paper Repeatability Package

This folder is the experimental repeatability package for the ATVA 2026
paper.

It bundles five end-to-end verification examples covering the major
new capabilities introduced by NNV3 — fairness, probabilistic, video,
graph, and weight perturbation — each runnable from a single command,
with all models and data shipped in-tree.

## Examples at a glance

| Folder | Capability | GPU? | Default sweep time |
|---|---|---|---|
| [GNNV](GNNV/) | Graph NN power-flow verification (PyG, GraphStar) | no | ~3.5 min |
| [FairNNV](FairNNV/) | Fairness on Adult-Income classifier | no | ~2 min |
| [ModelStar](ModelStar/) | Weight-perturbation on MNIST MLP | no | ~70 min |
| [VideoStar](VideoStar/) | 3D-CNN video verification (ZoomIn-4f, VolumeStar) | no | ~12 min |
| [ProbVer](ProbVer/) | Probabilistic verification of TinyYOLO (cp-star) | recommended | ~8 min (CPU) |

Total wall-clock for the full sweep on a single workstation:
**~95 minutes** (~1.5 hours).

## TL;DR — running the full paper reproduction

After the prerequisites and one-time setup below:

```bash
docker run --rm --gpus all \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv/code/nnv/examples/NNV3.0 \
    nnv3.0 ./run_all.sh --full
```

When it finishes, eyeball-compare your numbers against
[`EXPECTED_RESULTS.md`](EXPECTED_RESULTS.md) (which documents the
paper's headline numbers per example).

Read on for prerequisites, build, one-time setup, and the smaller
options (smoke run, single example, no-GPU).

---

## 1. Prerequisites

- **Docker** with the `--gpus` runtime if you want ProbVer
  GPU-accelerated.
- **MATLAB license**. The Dockerfile references a network license
  server via the `LICENSE_SERVER` build-arg (default Vandy's). External
  evaluators must override this — see *License setup* below.
- **NVIDIA GPU** with CUDA-capable drivers — recommended for ProbVer
  only (other examples are CPU). ProbVer falls back to CPU if no
  compatible GPU is detected.
- **Disk space**: the `nnv3.0` Docker image is ~26 GB (MATLAB R2024b
  dominates).

## 2. Build the Docker image

From the repo root:

```bash
docker build \
    --build-arg LICENSE_SERVER=27009@your.license.server \
    -t nnv3.0 \
    -f code/nnv/examples/NNV3.0/Dockerfile \
    .
```

The build:
- Installs MATLAB R2024b plus the toolboxes NNV needs.
- Copies the repo into the image at `/home/matlab/nnv`.
- Creates a Python venv and installs `requirement.txt` (used by ProbVer).
- Runs `install.m` to put NNV on the saved MATLAB path.

### License setup

Set `LICENSE_SERVER` to a license server you have access to:

```
--build-arg LICENSE_SERVER=27009@licenseserver.it.vanderbilt.edu   # internal default
--build-arg LICENSE_SERVER=27000@your-floating-server              # your institution
```

If you don't have a network license, you'll need to either use a
MathWorks-hosted MATLAB Online environment, or modify the Dockerfile
to use a different licensing strategy (e.g. mounting an
individual-license `network.lic`).

## 3. One-time host setup (when bind-mounting your repo)

The default `docker run` command below bind-mounts your local repo
into the container at `/home/matlab/nnv`, which **shadows** two things
the image built during `docker build`: the Python `.venv` and the
writable results / decompression directories. You need to recreate
them on the host once:

```bash
# (a) Build .venv on the host (required by ProbVer's cp-star Python helpers).
mkdir -p .venv && chmod 777 .venv
docker run --rm \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv \
    nnv3.0 bash -c "python3 -m venv .venv && \
                    .venv/bin/pip install --no-cache-dir -r requirement.txt"

# (b) Make per-example results / decompression directories writable from
#     inside the container (the matlab user inside has a different UID
#     than your host user).
chmod 777 \
    code/nnv/examples/NNV3.0/results \
    code/nnv/examples/NNV3.0/{FairNNV,ProbVer,VideoStar,ModelStar,GNNV}/results \
    code/nnv/examples/NNV3.0/ProbVer/yolo_2023/{onnx,vnnlib} \
    code/nnv/engine/nn/Prob_reach/Temp_files_mid_run \
    code/nnv/examples/NNV3.0/ModelStar/runtime
```

Skip step (a) if you don't plan to run ProbVer. The other examples
don't call into Python.

If you build the image without bind-mounting (i.e. run inside the
image's own `/home/matlab/nnv`), neither step is needed — but then any
local edits won't be picked up unless you rebuild.

## 4. Run

The orchestrator [`run_all.sh`](run_all.sh) handles all five examples
through a uniform interface. Per-example logs land in
`results/run_all_<timestamp>/<example>.log`.

### Full reproduction (~95 min)

```bash
docker run --rm --gpus all \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv/code/nnv/examples/NNV3.0 \
    nnv3.0 ./run_all.sh --full
```

### Smoke run (~5 min total — quick "does it all run?" check)

```bash
docker run --rm --gpus all \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv/code/nnv/examples/NNV3.0 \
    nnv3.0 ./run_all.sh --smoke
```

Smoke configurations per example:

- **GNNV**: 3 graphs × 1 ε × GCN only (~10 s)
- **FairNNV**: AC-1 only, 10 observations, ε ∈ {0.01} (~30 s)
- **ModelStar**: 1 layer (`fc_6`) × 4 fractions × 100 images (~1 min)
- **VideoStar**: 1 sample × ε=1/255 (~30 s)
- **ProbVer**: 1 instance, `nRand=50` (~3 min on CPU)

### No-GPU smoke (~2 min — skips ProbVer)

```bash
docker run --rm \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv/code/nnv/examples/NNV3.0 \
    nnv3.0 ./run_all.sh --smoke --no-gpu
```

### Single example

```bash
./run_all.sh --only gnnv          # gnnv | fairnnv | modelstar | videostar | probver
./run_all.sh --only gnnv --full   # full sweep for one example
```

### Direct (skip orchestrator)

Every example also runs standalone. From the NNV repo root:

```bash
matlab -batch "cd code/nnv/examples/NNV3.0/<Folder>; <entry-point>"
```

Per-example entry points:

| Folder | Entry point |
|---|---|
| GNNV | `run_gnn_experiments` |
| FairNNV | `run_fairnnv` |
| ModelStar | `run_expt_for_compute` |
| VideoStar | `run_zoomin_4f` |
| ProbVer | `run_probver` |

See each subfolder's `README.md` for full configuration options.

## 5. Outputs

Each example writes to its own `<Folder>/results/<timestamp>/`
subfolder. The orchestrator's per-run logs live under
`results/run_all_<timestamp>/<example>.log`. Outputs include CSVs,
`.mat` summaries, MATLAB diaries, and (FairNNV only) PNG/PDF figures
plus LaTeX tables.

To copy results out of a running container:

```bash
docker cp <container_id>:/home/matlab/nnv/code/nnv/examples/NNV3.0/<Folder>/results ./<Folder>_results
```

## 6. Verifying your reproduction matches the paper

After `--full` completes, eyeball-compare against
[`EXPECTED_RESULTS.md`](EXPECTED_RESULTS.md). It documents the paper's
headline numbers (Tables 3–6, Figure 5) per example, with the observed
verdicts and timings from a known-good reference run on this artifact.

For exact-match metrics (verified counts, SAT/UNSAT verdicts) the
artifact reproduces the paper's numbers cell-by-cell. Timings vary
with hardware.

## 8. Layout

```
NNV3.0/
├── README.md                   This file
├── run_all.sh                  Top-level orchestrator (--smoke|--full|--no-gpu|--only)
├── EXPECTED_RESULTS.md         Paper headline numbers per example
├── Dockerfile                  Build image (parameterized LICENSE_SERVER)
├── FairNNV/                    Fairness verification (Adult-Income)
├── GNNV/                       Graph NN verification (IEEE 24-bus)
├── ModelStar/                  Weight perturbation (MNIST MLP)
├── ProbVer/                    Probabilistic verification (TinyYOLO)
└── VideoStar/                  Video classification (ZoomIn-4f)
```
