# ToolComparison — NNV vs AIVL (ATVA 2026 artifact)

NNV vs MathWorks AI Verification Library (AIVL) head-to-head across the
six benchmarks reported in the ATVA 2026 paper. Smoke runs in ~5 min;
default runs in ~3 h on a reviewer-grade VM.

## System requirements

| | |
|---|---|
| MATLAB | R2025b (matches the CodeOcean capsule ceiling) |
| Toolboxes | Deep Learning Toolbox, Deep Learning Toolbox Converter for ONNX Model Format, Parallel Computing Toolbox |
| AIVL | MathWorks AI Verification Library Support Package (tarball recipe below) |
| RAM | ≥ 16 GB |
| Cores | ≥ 4 |
| Disk | ~250 MB (assets) + ~50 MB (results) |
| NNV | Cloned from https://github.com/verivital/nnv and `install.m` run once |

## One-command quickstart

```matlab
cd code/nnv/examples/NNV3.0/ToolComparison
run_toolcomparison('mode','smoke')      % ~5 minutes
run_toolcomparison                       % default mode, ~2.5 h
```

`run_toolcomparison` runs both halves (vnnlib + argmax) and renders
`tables/out/table_main.{tex,txt}`.

`run_toolcomparison` also honors the `TOOLCOMPARISON_MODE` env var
(`smoke|default|full`) so the NNV3.0 top-level `run_all.sh` orchestrator
can drive it without MATLAB-side arg plumbing. `'full'` is accepted as
an alias for `'default'`.

## Docker quickstart (R2025b)

```bash
cd code/nnv/examples/NNV3.0/ToolComparison
./build_image.sh                 # one-time, ~30–60 min
./run_all.sh                     # full sweep, all 6, nnv+aivl
./run_all.sh nnv                 # nnv only
./run_all.sh nnv,aivl rl oval21  # nnv+aivl, specific benchmarks
```

`build_image.sh` builds `nnv3.0:r2025b` with the Vanderbilt MATLAB
license server baked in. `run_all.sh` runs the sweep sequentially inside
the container, logging to `logs/` and writing per-benchmark results to
`results/`.

## NNV engine fixes included in this branch

Two small engine fixes are required for the grid to run cleanly on R2025b.
Both are committed directly on this branch (no separate install step) and
are intended to land in NNV `main`:

| File | Fix |
|---|---|
| `code/nnv/engine/nn/layers/FeatureInputLayer.m` | Reverts an Oct-2025 regression in `reach_star_single_input` that crashed on any `ImageStar` input (affects every FC-imported benchmark — acas_xu_p3/p4 and rl all hit this). Restores the pre-Oct-2025 `affineMap([], -obj.Mean)` form that works for both `Star` and `ImageStar`. |
| `code/nnv/engine/nn/NN.m` | `start_pool` bails when `getCurrentTask()` indicates a parfeval worker context, so `exact-star`'s pool isn't re-created from inside a worker. Also switches `parpool('local',n)` to a `parcluster()` handle to bypass the default 8-worker Processes-profile cap. Required for any `exact-star` call wrapped in `parfeval` (acas_xu_p3/p4 + rl). |

If you're running this artifact from a fresh clone of NNV `main` that
*doesn't* yet have these fixes merged, `cd code/nnv && git log
code/nnv/engine/nn/layers/FeatureInputLayer.m code/nnv/engine/nn/NN.m`
will tell you whether the fixes are present. If absent, cherry-pick this
branch's engine commit or apply the diffs manually.

## Install AIVL

```matlab
cd code/nnv/examples/NNV3.0/ToolComparison/utils
aivl_install
```

This finds the tarball at `../atva26-aivl.tar.gz` (symlinked to
`utils/atva26-aivl.tar.gz`), extracts the AIVL Support Package into
`userpath()/SupportPackages/`, and persists the addpath via `startup.m`.
Restart MATLAB once after the first install, then verify:

```matlab
which verifyNetworkRobustness   % should resolve to the aivnv dir
which estimateNetworkOutputBounds
```

**No AIVL license?** Run NNV-only:

```matlab
run_toolcomparison('mode','default','tools',{'nnv'})
```

The harness works unchanged; AIVL rows are simply absent from `table_main`.

## Benchmark catalog

| Benchmark | N inst | Network | Property | NNV (default) | AIVL |
|---|---:|---|---|---|---|
| **acas_xu_p3** | 20 | ACAS Xu FC ReLU (6 layers) | half-space VNNLIB (prop_3) | `approx-star`, `exact-star`, `relax-star-range-50` | `estimate-bounds` |
| **acas_xu_p4** | 20 | ACAS Xu FC ReLU (6 layers) | half-space VNNLIB (prop_4) | `approx-star`, `exact-star`, `relax-star-range-50` | `estimate-bounds` |
| **rl** | 50 | Cartpole + Lunarlander FC (tanh+ReLU) | half-space VNNLIB | `approx-star`, `exact-star`, `relax-star-range-50` | `estimate-bounds` |
| **oval21** | 30 | CIFAR ResNet (3 nets × 10 props) | half-space VNNLIB | `approx-star` | `estimate-bounds` |
| **collins_rul** | 62 | Small 1D Conv (Collins RUL) | half-space VNNLIB | `approx-star` | `estimate-bounds` |
| **mnist_resnet8** | 25 imgs × 4 ε | Custom MNIST ResNet-8 (additionLayer) | argmax robustness (L∞) | `relax-star-area-50` | `deep-poly` (verifyNetworkRobustness) |

In **smoke mode** every benchmark runs **one instance with one NNV
algorithm only**; AIVL is skipped (avoids the ~30 s AIVL warmup × 6
benchmarks). Smoke produces 6 result rows in ~5 minutes.

In **default mode** every benchmark runs every algorithm above on its
full instance count. AIVL rows are included.

### Notes on algorithm choices

- **FC nets (ACAS, RL):** `approx-star` + `exact-star` + `relax-star-range-50`. The full FC grid matches NNV's CAV'23 convention; exact-star is feasible on these small networks.
- **CNN / ResNet (oval21, collins_rul, mnist_resnet8):** `approx-star` for oval21/collins_rul; `relax-star-area-50` for mnist_resnet8. `exact-star` never finishes under VNN-COMP timeouts on Conv nets. `relax-star-area-50` produced identical V/X counts to `approx-star` on prior oval21/collins_rul runs, so it's omitted for those two as redundant. mnist_resnet8 uses `relax-star-area-50` per its argmax robustness convention.
- **Other `relax-star-{range,area}-25/-75/-100` variants:** deliberately omitted. Prior data shows `-25` ≈ `approx-star` (no value-add) and `-100` ≈ pure interval (loses precision). `-50` is the documented sweet spot.
- **AIVL `deep-poly`** is the right pick for `mnist_resnet8` (argmax-form, test set carries `ytrue`). The five VNNLIB-style benchmarks use `estimate-bounds` instead — VNNLIB-form `verifyNetworkRobustness` is an R2026a feature.

## Mode reference

| Mode | What runs | Wall-clock | Disk written |
|---|---|---|---|
| `'smoke'` | 1 inst × 1 NNV alg per benchmark, AIVL skipped | ~5 min | <1 MB results |
| `'default'` | Full grid above, full instance counts | ~2.5 h | ~5 MB results |
| `'full'` | Alias for `'default'` (paper convention) | ~2.5 h | ~5 MB results |

## Reading the results

| Artifact | Where |
|---|---|
| Per-benchmark raw rows | `results/<bench>.mat` — variable `results`, MATLAB table, schema in `utils/tool_utils.m` |
| Consolidated table | `tables/out/table_main.tex` (LaTeX) + `table_main.txt` (plaintext) |
| Per-instance log | console only (resume-skip mechanism uses `tool_utils.has_instance`) |
| Per-benchmark Docker log | `logs/<bench>.log` |
| Run-wide issue summary | `ISSUES.md` (auto-appended by `run_all.sh`) |

Table columns: V (verified) / X (violated, i.e. counterexample) / ? (unknown) / T/O (timeout) / Err / Mean t / PAR-2.

PAR-2 is the VNN-COMP scoring metric — unsolved instances contribute 2 × their timeout to the mean.

## Reproducibility notes

- Reruns are **resume-safe**. `(tool, benchmark, instance_id, algorithm)` quadruples already present in the `.mat` are skipped.
- To force a re-run of timed-out rows under a higher timeout, purge them first:
  ```matlab
  u = tool_utils;
  u.purge_status('results/acas_xu_p3.mat', 'timeout', 'nnv', 'exact-star');
  run_toolcomparison;
  ```
- The harness creates and destroys the parpool between instances on the
  vnnlib half (NNV-AIVL worker state isolation); the argmax half runs
  inline (mnist_resnet8 doesn't need cross-tool isolation and the inline
  path is ~6× faster).
- The smoke `.mat` files are a strict subset of the default ones —
  running smoke then default is supported and just adds rows.

## Known limitations / non-goals

- No α/β-CROWN AIVL bridge — that lands in R2026a, out of scope for the
  R2025b ceiling.
- VNNLIB-form `verifyNetworkRobustness` (true `deep-poly` on half-space
  output specs) is also R2026a-only. The five VNNLIB benchmarks use
  `estimate-bounds` instead, which is the documented R2025b AIVL
  VNNLIB path.
- No CIFAR2020 — NNV's existing data shows it loses badly to AIVL on
  CIFAR2020; out of scope for the head-to-head story.
- No cgan / TransposedConv2D — AIVL R2025b doesn't support
  `TransposedConv2D` or custom `FunctionLayer`, blocking the cgan
  comparison entirely.
- No dubinsrejoin from the RL bench — R2025b ScalingLayer folding
  produces a Star/affineMap dim mismatch we couldn't fix without
  engine-level changes.

## Repository layout

```
ToolComparison/
├── README.md                          ← this file
├── run_toolcomparison.m               ← MATLAB entry (smoke|default|full)
├── run_all.sh                         ← Docker-driven full sweep
├── run_step_by_step.sh                ← interactive per-benchmark runner
├── build_image.sh                     ← build nnv3.0:r2025b
├── drivers/
│   ├── run_vnnlib_half.m              5 VNNLIB benchmarks
│   └── run_argmax_half.m              mnist_resnet8
├── benchmarks/                        bundled ONNX + VNNLIB assets
│   ├── acas_xu_p3/, acas_xu_p4/, rl/, oval21/, collins_rul/, mnist_resnet8/
├── utils/
│   ├── tool_utils.m                   schema, append, has_instance, par2
│   ├── reach_opt_for.m                algorithm string → reachOpt
│   ├── rebuild_for_aivl.m             strip ScalingLayer/VerifyBatchSize
│   ├── aivl_install.m                 AIVL tarball installer
│   ├── toolbox_install.m              tarball extractor
│   ├── parse_argmax_vnnlib.m
│   ├── load_mw_network.m
│   ├── addpath_shared.m
│   └── atva26-aivl.tar.gz             AIVL Support Package tarball
├── tables/
│   ├── make_table_main.m              consolidated table renderer
│   └── out/                           generated table_main.{tex,txt}
├── results/                           created at runtime
├── logs/                              created by run_all.sh
└── ISSUES.md                          auto-appended summary
```

## Citation

```bibtex
@inproceedings{tumlin2026nnv3,
  title  = {NNV 3.0 Tool Paper},
  author = {Tumlin, Anne M. and Manzanas Lopez, Diego and Johnson, Taylor T. and others},
  booktitle = {ATVA 2026},
  year   = {2026}
}
```
