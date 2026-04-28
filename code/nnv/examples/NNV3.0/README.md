# NNV3 — Tool Paper Repeatability Package

Experimental repeatability package for the ATVA 2026 NNV3 paper.
Five end-to-end verification examples covering NNV3's major new
capabilities: fairness, probabilistic, video, graph, and weight
perturbation.

> **Cite**: Tumlin, Sasaki, et al. *NNV3: Expanding Neural Network
> Verification to New Architectures and Domains.* ATVA 2026.
>
> **Issues / questions**: https://github.com/atumlin/nnv/issues

| Folder | Capability | GPU? | Smoke | Default sweep |
|---|---|---|---:|---:|
| [GNNV](GNNV/) | Graph NN power-flow verification (PyG, GraphStar) | no | ~10 s | ~3.5 min |
| [FairNNV](FairNNV/) | Fairness on Adult-Income classifier | no | ~30 s | ~2 min |
| [ModelStar](ModelStar/) | Weight-perturbation on MNIST MLP | no | ~1 min | ~70 min |
| [VideoStar](VideoStar/) | 3D-CNN video verification (ZoomIn-4f, VolumeStar) | no | ~30 s | ~12 min |
| [ProbVer](ProbVer/) | Probabilistic verification of TinyYOLO (cp-star) | recommended | ~3 min | ~8 min (CPU) |

**Total full-sweep wall time** on a single workstation: **~95 minutes**
(~1.5 hours). Smoke-only across all five: **~5 minutes**.

```
NNV3.0/
├── README.md                  This file
├── EXPECTED_RESULTS.md        Paper headline numbers vs. observed (per example)
├── Dockerfile                 License-agnostic build recipe
├── run_all.sh                 Top-level orchestrator
├── FairNNV/                   ↘
├── GNNV/                       │  Each subfolder has its own README.md
├── ModelStar/                  │  with overview, references, configuration,
├── ProbVer/                    │  outputs, and observed runtime.
└── VideoStar/                 ↗
```

---

## Quickstart (~30 minutes to a working smoke run)

Five steps, in order. Each builds on the previous. By the end you'll
have run a smoke configuration of every example and verified the
artifact works on your machine. After that, the *Full reproduction*
section kicks off the paper-faithful sweep.

### Step 1 — Prerequisites

You need:

- **Docker** (any recent version) with the `--gpus` runtime if you
  want ProbVer GPU-accelerated.
- **MATLAB R2024b license** at run time (not at build time). Three
  reviewer patterns are supported in Step 3 below; pick one.
- **NVIDIA GPU** with CUDA-capable drivers — recommended for ProbVer
  only (other examples are CPU). ProbVer falls back to CPU if no
  compatible GPU is detected.
- **Disk space**: ~26 GB for the Docker image; another few GB for
  per-run outputs.

Confirm Docker is working:

```bash
docker --version                         # any recent version
docker info | grep -i runtime            # should list 'nvidia' if you want GPU
```

If you'll use a network license server (Pattern A in Step 3), confirm
it's reachable from your build host before going further:

```bash
nc -zv your.license.server 27000         # substitute your port@host
```

### Step 2 — Build the Docker image

The build is **license-agnostic** — MATLAB is installed but never
launched, so no license check happens at build time. From the **NNV
repo root** (the directory containing `code/`):

```bash
docker build -t nnv3.0 -f code/nnv/examples/NNV3.0/Dockerfile .
```

What the build does:

- Installs MATLAB R2024b plus the toolboxes NNV needs (~15-20 min the first time).
- Copies the repo into the image at `/home/matlab/nnv`.
- Creates a Python venv inside the image and installs `requirement.txt` (used by ProbVer).

NNV's MATLAB-path setup is **not** baked into the image — each example
self-adds NNV at run time, which is why the build never needs to launch
MATLAB.

Verify:

```bash
docker images | grep nnv3.0              # should show nnv3.0:latest, ~24 GB
```

### Step 3 — License setup (pick one)

This artifact supports three reviewer scenarios. The example commands
in the rest of this README assume **Pattern A** (network license
server). If you're using Pattern B, swap the `-e MLM_LICENSE_FILE=...`
flag for `-v /path/to/network.lic:/opt/matlab/R2024b/licenses/network.lic:ro`.

#### Pattern A — Institutional network license server

Most common for academic / corporate MATLAB users. Pass `port@hostname`
as an env var at every `docker run`:

```bash
docker run --rm \
    -e MLM_LICENSE_FILE=27000@your.license.server \
    nnv3.0 matlab -nodisplay -batch "disp('license OK'); exit"
```

(Replace `27000@your.license.server` with what your institution
gave you — format is `port@hostname`.)

#### Pattern B — Individual MathWorks license file

Mount your `network.lic` (from your activated MathWorks account or a
floating-license file) into the container's licenses directory:

```bash
docker run --rm \
    -v /path/to/network.lic:/opt/matlab/R2024b/licenses/network.lic:ro \
    nnv3.0 matlab -nodisplay -batch "disp('license OK'); exit"
```

#### Pattern C — No local MATLAB license

Three options:

- **CodeOcean capsule** (recommended for hosted evaluation): the paper
  references a [hosted evaluation capsule](https://codeocean.com/capsule/1928763/tree)
  that runs MATLAB in CodeOcean's environment with no local license.
- **MathWorks 30-day trial**: register for a [free trial license](https://www.mathworks.com/campaigns/products/trials.html)
  and use Pattern B above with the trial's `network.lic`.
- **Read-only inspection**: this folder ships
  [`EXPECTED_RESULTS.md`](EXPECTED_RESULTS.md) with the verdicts and
  timings observed when the artifact was last validated end-to-end.
  Reviewers without MATLAB can eyeball-compare those tables against
  the paper's Tables 3-6 / Figure 5 directly without running anything.

### Step 4 — One-time host setup (when bind-mounting your repo)

The `docker run` commands below bind-mount your local repo into the
container at `/home/matlab/nnv`. The bind-mount **shadows** two things
the image built during `docker build`: the Python `.venv` and writable
results / decompression directories. Recreate them on the host once,
from the **NNV repo root**.

> The `chmod 777` block uses bash brace expansion. If you're in a
> `dash`/POSIX `sh`, run `bash` first or expand the paths manually.

```bash
# (a) Build .venv on the host (required only by ProbVer's Python helpers).
mkdir -p .venv && chmod 777 .venv
docker run --rm \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv \
    nnv3.0 bash -c "python3 -m venv .venv && \
                    .venv/bin/pip install --no-cache-dir -r requirement.txt"

# (b) Make per-example results / decompression directories writable from
#     inside the container. The container's matlab user (UID 1001) has
#     a different UID than your host user, so it can't write into
#     host-owned directories without these chmods.
chmod 777 code/nnv/examples/NNV3.0/results
chmod 777 code/nnv/examples/NNV3.0/FairNNV/results
chmod 777 code/nnv/examples/NNV3.0/GNNV/results
chmod 777 code/nnv/examples/NNV3.0/ModelStar/results
chmod 777 code/nnv/examples/NNV3.0/ModelStar/runtime
chmod 777 code/nnv/examples/NNV3.0/ProbVer/results
chmod 777 code/nnv/examples/NNV3.0/ProbVer/yolo_2023/onnx
chmod 777 code/nnv/examples/NNV3.0/ProbVer/yolo_2023/vnnlib
chmod 777 code/nnv/examples/NNV3.0/VideoStar/results
chmod 777 code/nnv/engine/nn/Prob_reach/Temp_files_mid_run
```

Skip step (a) if you don't plan to run ProbVer; the other examples
don't call into Python.

If you don't bind-mount (i.e. run inside the image's own
`/home/matlab/nnv`), neither step is needed — but then any local
edits won't be picked up unless you rebuild.

### Step 5 — Run the smoke (~5 minutes)

Quick end-to-end check of all five examples on small per-example
configurations. From the **NNV repo root**:

```bash
docker run --rm --gpus all \
    -e MLM_LICENSE_FILE=27000@your.license.server \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv/code/nnv/examples/NNV3.0 \
    nnv3.0 ./run_all.sh --smoke
```

Variants:

- **No GPU**: drop `--gpus all` and add `--no-gpu` to skip ProbVer.
- **Subset of examples**: `--only gnnv` (or `fairnnv|modelstar|videostar|probver`).
- **Full reproduction** instead of smoke: swap `--smoke` for `--full` (~95 min, see next section).

What you should see (sample, abbreviated):

```
================================================================
  NNV3.0 — run_all  (mode=smoke, no_gpu=0, only=all)
================================================================
[GNNV] >> matlab -batch "..."
   OK
[FairNNV] >> matlab -batch "..."
   OK
[ProbVer] >> matlab -batch "..."
   OK
[VideoStar] >> matlab -batch "..."
   OK
[ModelStar] >> matlab -batch "..."
   OK
================================================================
  All requested examples completed.
================================================================
```

If any example FAILs, the orchestrator continues with the rest and
reports which example failed at the end. Inspect the failing
example's log under `results/run_all_<timestamp>/<example>.log`; see
*Troubleshooting* below for common issues.

---

## Full reproduction (~95 minutes)

Same command pattern as Step 5, with `--smoke` swapped for `--full`:

```bash
docker run --rm --gpus all \
    -e MLM_LICENSE_FILE=27000@your.license.server \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv/code/nnv/examples/NNV3.0 \
    nnv3.0 ./run_all.sh --full
```

ModelStar dominates the wall time (~70 min of the ~95 min total) —
`fc_4` perturbations alone take ~30 minutes because deeper perturbed
layers compose with more downstream ReLU operations. See the
at-a-glance table at the top for per-example timing.

---

## Running a single example

Most reviewers will use:

```bash
./run_all.sh --only gnnv          # one of: gnnv | fairnnv | modelstar | videostar | probver
./run_all.sh --only modelstar --full
```

Per-example entry points (if you ever want to skip the orchestrator
and call MATLAB directly):

| Folder | Entry point |
|---|---|
| GNNV | `run_gnn_experiments` |
| FairNNV | `run_fairnnv` |
| ModelStar | `run_expt_for_compute` |
| VideoStar | `run_zoomin_4f` |
| ProbVer | `run_probver` |

See each subfolder's `README.md` for full configuration options
(custom epsilons, subset of architectures, etc.).

---

## Outputs and verification

### Where outputs land

- Each example writes to `<Folder>/results/<timestamp>/` containing
  CSVs, `.mat` summaries, MATLAB diaries, and (FairNNV only) PNG/PDF
  figures plus LaTeX tables.
- The orchestrator writes per-example MATLAB stdout/stderr to
  `results/run_all_<timestamp>/<example>.log`.

### Comparing against the paper

After a `--full` run, eyeball-compare against
[`EXPECTED_RESULTS.md`](EXPECTED_RESULTS.md). It contains side-by-side
tables of paper numbers vs. observed numbers (from a known-good
validation run of this artifact) for every example.

For verdicts (verified counts, SAT/UNSAT, fair/unfair %), the
artifact reproduces the paper's headline numbers cell-by-cell.
Timings vary with hardware — expect ±50% relative drift.

### Copying outputs out of a stopped container

Only relevant if you ran without `-v $PWD:...` bind-mount. The
recommended setup uses bind-mount, so outputs already appear directly
on your host under `code/nnv/examples/NNV3.0/<Folder>/results/`.

```bash
docker cp <container_id>:/home/matlab/nnv/code/nnv/examples/NNV3.0/<Folder>/results ./<Folder>_results
```

---

## Troubleshooting

### Common (most reviewers will hit at least one)

| Symptom | Cause | Fix |
|---|---|---|
| `Sign-in failed` / MATLAB asks for credentials | No license configured at `docker run` time | Add `-e MLM_LICENSE_FILE=...` (Pattern A) or mount a license file (Pattern B). See *Step 3*. |
| `mkdir: Permission denied` on a `results/` folder | Container's `matlab` user (UID 1001) ≠ host user UID, can't write | Run *Step 4 (b)* — the `chmod 777` block — once on the host |
| `Python virtual environment not found at /home/matlab/nnv/.venv` | Bind-mount has shadowed the image's pre-built venv | Run *Step 4 (a)* — the venv-build block |

### GPU / Python edge cases

| Symptom | Cause | Fix |
|---|---|---|
| ProbVer falls back to CPU even with `--gpus all` | torch wheels in `.venv` require CUDA 13; older drivers support only CUDA 12 | Acceptable (slower but works). For GPU acceleration, rebuild the venv with `pip install torch --index-url https://download.pytorch.org/whl/cu121` |
| GNNV `Unrecognized function or variable 'glpk'` | Bind-mount shadows the image's tbxmanager/glpkmex; linprog falls back to glpk at the loosest ε | Already handled gracefully for GNNV (degrades to "unknown"). Elsewhere: rebuild the image without bind-mounting. |
| `Unable to write file ... Direction_data.mat` | UID mismatch (same root cause as `mkdir: Permission denied`) but on ProbVer's intermediate-files dir | Same fix — *Step 4 (b)* covers `Prob_reach/Temp_files_mid_run` |

### Behavioral

| Symptom | Cause | Fix |
|---|---|---|
| One example fails — does the orchestrator stop? | No — `run_all.sh` wraps each example so failure doesn't abort the rest | Per-example exit status is printed; orchestrator returns non-zero if any failed |
| Stale `.onnx` after rebuilding the model | `run_probver.m` only decompresses if the `.onnx` is missing | Delete `ProbVer/yolo_2023/onnx/TinyYOLO.onnx` to force fresh decompression |
| First build is slow | First-time `mpm install` downloads ~10 GB of MATLAB toolboxes | Normal (~15-20 min). Layer cache makes subsequent builds nearly instant. |
