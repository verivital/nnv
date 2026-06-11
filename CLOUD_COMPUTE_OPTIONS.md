# Cloud / Parallel Compute Options for NNV — Research Report

_Generated 2026-06-08 via deep-research (24 sources, 116 claims extracted, 25 adversarially verified, 19 confirmed). Goal: fan out the embarrassingly-parallel NNV test suite (~833 tests, ~70 min) and VNN-COMP 2025 sweep (dozens of independent instances) to get results in minutes, not hours. CPU/LP-bound, no GPU._

## TL;DR
- **NNV is already a PUBLIC repo (`github.com/verivital/nnv`) → the fastest, free, zero-licensing path is GitHub Actions matrix builds with `matlab-actions/setup-matlab` + `run-tests`.** On **public** repos these actions **auto-license ALL products for free** (except transformation products like Coder/Compiler, which NNV doesn't use) — so all four NNV toolboxes (Deep Learning, Parallel Computing, Optimization, Statistics & ML) work with **no token and no campus-license config**. A `strategy.matrix` shards work across parallel jobs (≤256 jobs/run, ~20 concurrent runners on free/Team tiers). **Hours → minutes.**
- For heavy/memory-bound exact-star or very large sweeps: **ACCRE SLURM array jobs** reusing Vanderbilt's campus license server (no batch token, no admin status needed) — strong but **needs ACCRE-specific confirmation** (low-confidence/inferred below).
- **Docker on cloud batch (AWS/GCP/Azure)** is the most scalable but has the **worst licensing friction** — avoid as a first step.

## Ranked options

### 1. GitHub Actions matrix + matlab-actions  ⭐ PRIMARY (free, lowest setup)
- **How:** `.github/workflows/*.yml` with `uses: matlab-actions/setup-matlab@v2` then `uses: matlab-actions/run-tests@v2`, plus `strategy: matrix:` to shard tests/instances into separate parallel jobs (Linux/Windows/macOS).
- **Licensing:** **PUBLIC repo → all products auto-licensed, free, no token.** (verbatim: "Public project — The … actions automatically license all products for you, except for transformation products". confidence: high, 3‑0)
- **Parallelism ceiling:** ≤256 jobs/run; ~20 concurrent runners on free/Team (paid/larger runners raise this). `use-parallel: true` adds intra-job PCT parallelism (also auto-licensed on public).
- **Cost:** Free for public repos on GitHub-hosted runners (within Actions minutes).
- **Gotchas:** NNV has **no `.prj`**, so `run-tests` falls back to "all tests in the repo root and subfolders" (still works, just not the `.prj` "Test label" default). You author a small dispatcher so each matrix shard runs only its slice.

### 2. ACCRE SLURM array jobs (free, reuses campus license)
- **How:** `sbatch --array=1-N`, each task runs `matlab -batch "run_instance(N)"` headlessly on a cluster node, checking out a concurrent seat from Vanderbilt's network license server over TCP 27000–27010.
- **Licensing:** Reuses the existing campus network license — **no batch token, no admin status, no cloud-licensing** needed.
- **Cost:** Free (campus allocation).
- **Confidence: LOW / inferred** — the general NLM mechanism is verified, but ACCRE's specific MATLAB module name, available **Optimization Toolbox / PCT concurrent seats**, and array-size/seat-contention limits were **not** verified. Confirm with ACCRE / the campus MATLAB admin. (ref: help.accre.vanderbilt.edu MATLAB + Job-Arrays pages.)

### 3. Private-repo CI: batch token OR self-hosted Linux runner
- If NNV CI must run on a **private** repo, matlab-actions auto-licenses **nothing**. Two paths:
  - **MathWorks batch licensing token** → set as repo secret, mapped to `MLM_LICENSE_TOKEN`. **But tokens are a PILOT** — eligibility form (name/email/license #), **not self-service**, and **academic/TAH eligibility is NOT publicly confirmed**; the form was intermittently "currently unavailable." Don't assume the Vanderbilt license grants one. (high, 3‑0)
  - **Self-hosted runner** with a traditionally-licensed MATLAB pointed at the campus license server. Must be **Linux/macOS** (setup-matlab does **not** support Windows self-hosted runners).

### 4. Docker `mathworks/matlab` on AWS/GCP/Azure Batch (most scalable, most friction)
- **Image gotcha:** official `mathworks/matlab` ships **bare core MATLAB, no toolboxes** (Ubuntu 24.04) — all four NNV toolboxes must be `mpm`-installed. Even `mathworks/matlab-deep-learning` bundles DL/PCT/Statistics but **lacks Optimization Toolbox** (NNV's `linprog`), so a **custom image build is unavoidable**. (high, 3‑0)
- **Licensing:** needs **either** a batch token **or** `-e MLM_LICENSE_FILE=27000@server` pointing at a network license manager. Standing up your **own** AWS network-license-manager (official CloudFormation ref-arch exists) requires **"Administrator status for the network license"** — a likely blocker (you're probably not the campus license admin). (high, 3‑0)
- **Verdict:** highest ceiling, but only worth it once #1/#2 are exhausted.

### 5. MATLAB Parallel Server / Cloud Center / MATLAB Online
- Parallel Server (incl. campus on-demand/AWS) scales `parfor`/`batch` to many workers but is a heavier setup and licensing-gated. MATLAB Online has **batch/automation limits** (interactive-oriented). Useful if Vanderbilt already runs Parallel Server, otherwise not the fast path.

### 6. CodeOcean — NOT a primary path
- You can **run** existing capsules via API but **cannot create capsules** via API (web-UI only). So fan-out would require first hand-building a parameterized MATLAB capsule. Deprioritized vs the public-repo path. (from project memory)

## Recommended PRIMARY + quick-start
**Make/refresh a matrix-sharded `run-tests` workflow on the public NNV repo** (zero licensing setup, free, fastest to minutes):
1. `.github/workflows/test-matrix.yml`:
   - `strategy: matrix: shard: [1,2,3,4,5,6,7,8]`
   - `uses: matlab-actions/setup-matlab@v2` (with `release: R2026a`)
   - `uses: matlab-actions/run-command@v2` calling a dispatcher `run_shard(${{ matrix.shard }}, 8)` that selects this shard's tests/VNN-COMP instances and runs them headlessly.
2. Same pattern for the VNN-COMP sweep: one matrix entry per benchmark category (or per instance), each a `matlab -batch` of `run_vnncomp_instance`.
3. Collect per-shard results as workflow artifacts; merge into one CSV/MD.

**If it must be private later:** apply for a batch token (and/or stand up a Linux self-hosted runner on the campus network). **For the heaviest exact-star runs:** ACCRE SLURM arrays.

## ✅ As implemented in this repo (matrix CI + high-memory fallback)

Added `.github/workflows/test-matrix.yml` + `code/nnv/tests/run_shard.m`:

- **Standard fan-out:** a `prepare` job builds a shard list from `num_shards` (default 8); the `standard` job runs `strategy.matrix.shard` on free `ubuntu-latest` runners, each calling `run_shard(k, N, 'standard')`. `run_shard` enumerates the whole suite (`TestSuite.fromFolder`), **round-robin partitions** it into shard k of N, runs it headlessly, writes a JUnit XML artifact, and `assertSuccess`-fails the job on any failure. Validated on R2026a: 1095 test elements → ~1033 standard (62 heavy excluded). 8 shards ≈ **~10 min wall-clock vs ~70 min sequential.**
- **Memory crux:** the free public-repo `ubuntu-latest` runner is **~16 GB RAM / 4 vCPU**. NNV tests that exceed this are matched by `highMemPatterns` in `run_shard.m` and **excluded from the standard shards** so CI stays green/fast (no OOM). Current heavy set = **62 elements / 13 files**: VolumeStar, Conv3D, AveragePooling3D, Segmentation/SEGNET/PixelClassification, exact-star NNCS, comprehensive_network, and `_gpu_` tests (no GPU on runners). **Edit `highMemPatterns` as the full-suite duration/OOM report reveals more.**

### High-memory fallback (ranked) — run the excluded tests here
Use `run_shard(1, 1, 'highmem')` (or shard it) on a bigger machine:

1. **Self-hosted runner (recommended primary fallback) — free, integrated.** Register a big-RAM machine (your R2026a workstation, or a lab Linux box with 64+ GB) as a GitHub self-hosted runner labeled `nnv-highmem`, set repo **variable `NNV_HIGHMEM_RUNNER=nnv-highmem`**, then run the workflow with **`run_highmem=true`** — the `highmem` job routes there automatically and reuses your already-activated campus license.
   - ⚠️ `setup-matlab` cannot *install* MATLAB on **Windows** self-hosted runners. On Windows, use a **pre-installed** R2026a (skip `setup-matlab`; the `run-command` step uses the local MATLAB), or use a **Linux/macOS** self-hosted runner.
2. **ACCRE SLURM array (best for the heaviest / many large jobs) — free, scales.** `sbatch --array=1-N` with each task `matlab -batch "cd code/nnv; install; addpath(genpath('tests')); run_shard($SLURM_ARRAY_TASK_ID, N, 'highmem')"` on high-memory nodes, checking out concurrent seats from the campus license server. (Confirm ACCRE MATLAB module + Optimization/PCT seats — see open questions.)
3. **GitHub larger runners (paid, simplest).** If you have Team/Enterprise, set `NNV_HIGHMEM_RUNNER` to a larger-runner label (32–64 GB) — no infra to manage, but billed per-minute.

**Where it runs:** the workflow triggers on `workflow_dispatch` (manual, with `release`/`num_shards`/`run_highmem` inputs), `pull_request` to `master` (so PR #290 exercises it on `verivital/nnv`), and `push` to `transformer`. To run it on your fork, **enable Actions on `ttj/nnv`** (forks have Actions off by default) or trigger via the Actions tab.

### Empirical CI results (first runs, 2026-06-08, ttj/nnv)

- **It works.** The matrix executes correctly: `prepare` → 8 parallel shards → tests → JUnit + `assertSuccess`. Public-repo auto-licensing worked on the fork (no token).
- **Timing / speedup:** local serial = **33.9 min** (1095 tests). Matrix (8 shards): **~19 min cold** (first run, full `setup-matlab` install per shard), **~13.5 min on the next run** (install cached) — **~2.5×**, *including* the heavy tests. Bottlenecks: (a) ~5-min cold MATLAB install per shard (cached after the first run), and (b) load imbalance (one shard drew the slow tests). Levers to reach ~5–6×: install caching (automatic on re-runs), **duration-balanced sharding** (bin-pack by the per-test durations in `full_suite_results_r2026a.csv`), and more shards.
- **High-memory reality (key finding):** of 60 anticipated-heavy elements, **45 ran fine on the 16 GB runner**; only the few in the 2 crashed shards actually exceed it (MATLAB **exit 255** = fatal crash, *not* a clean 137 OOM-kill). So "anticipated-heavy" was mostly wrong — **the vast majority need no special handling.** Only a tiny number go to the high-mem fallback.
- **Diagnosing crashes:** GitHub discards a crashed job's stdout, so `run_shard` registers `ProgressPlugin`, which appends a *flushed* `START`/`DONE` line per test; on a crashed shard the last `START` without a `DONE` names the exact culprit (uploaded as an artifact). Reusable for any future CI crash.
- **Standard-shard exclusions (minimal):** GPU tests (no GPU on runners) + 2 tests broken under `runtests` (SLM `%%`-scoping) + the confirmed crashers only. Everything else — incl. ~58 of 62 "heavy" — runs in standard CI.
- **Other CI-only failures:** 4 `test_all_tutorial` sections (NNCS AEBS / InvertedPendulum, MNIST verify, SetRepresentations) fail on the Linux runner but **pass locally** — a CI-environment difference (data/products/platform), not a code bug; to investigate.

## Open questions to confirm
1. Does Vanderbilt's TAH/Campus-Wide license qualify for a MATLAB **batch token**, and is the pilot form open to academics? (unconfirmed)
2. ACCRE: is there a MATLAB module reaching the campus license server from compute nodes, and how many **Optimization Toolbox / PCT concurrent seats** + array limits? (unverified)
3. For a private repo, can a Linux self-hosted runner reach the campus license server (IT policy / network)?
4. Can you (or a delegate) get **network-license administrator** status + "configured for cloud use" (gates the AWS NLM / cloud-container concurrent-license paths)?

## Caveats
- A general MathWorks "licensing in the cloud" page claimed TAH licenses are cloud-ready; those specific claims were **REFUTED** in verification (0‑0 / 1‑0). The **Cloud Center container doc** is authoritative for the container case; treat TAH cloud eligibility as **needs admin confirmation**.
- Batch-token availability is **time-sensitive** (pilot).

## Key sources
- matlab-actions/setup-matlab + run-tests (primary): github.com/matlab-actions/setup-matlab , /run-tests
- MathWorks batch tokens (primary): mathworks.com/support/batch-tokens.html
- MATLAB container on Docker Hub / Cloud Center (primary): mathworks.com/help/cloudcenter/ug/matlab-container-on-docker-hub.html
- License manager for MATLAB on AWS (primary): github.com/mathworks-ref-arch/license-manager-for-matlab-on-aws
- ACCRE MATLAB + Job Arrays (primary): help.accre.vanderbilt.edu (MATLAB_on_the_ACCRE_Cluster ; Parallel_Processing_and_Job_Arrays)
