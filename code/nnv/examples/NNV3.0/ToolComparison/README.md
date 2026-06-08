# ToolComparison — NNV vs AIVL (ATVA 2026 artifact)

NNV vs MathWorks AI Verification Library (AIVL) head-to-head across the
six benchmarks reported in the ATVA 2026 paper. Smoke runs in ~1 min
(parpool-startup-dominated); default runs in ~3 h on a reviewer-grade VM.

## System requirements

| | |
|---|---|
| MATLAB | R2025b |
| Toolboxes | Deep Learning Toolbox, Deep Learning Toolbox Converter for ONNX Model Format, Parallel Computing Toolbox |
| AIVL | MathWorks AI Verification Library Support Package (auto-installed via `mpm` in both Docker flows; see "Install AIVL" below for host-MATLAB users) |
| RAM | ≥ 16 GB |
| Cores | ≥ 4 |
| Disk | ~250 MB (assets) + ~50 MB (results) |
| NNV | Cloned from https://github.com/verivital/nnv and `install.m` run once |

## One-command quickstart

```matlab
cd code/nnv/examples/NNV3.0/ToolComparison
run_toolcomparison('mode','smoke')      % ~1 minute
run_toolcomparison                       % default mode, ~2.5 h
```

`run_toolcomparison` runs both halves (vnnlib + argmax) and renders
`tables/out/table_main.{tex,txt}`.

`run_toolcomparison` also honors the `TOOLCOMPARISON_MODE` env var
(`smoke|default|full`) so the NNV3.0 top-level `run_all.sh` orchestrator
can drive it without MATLAB-side arg plumbing. `'full'` is accepted as
an alias for `'default'`.

## Docker quickstart

The artifact-evaluation Docker path is documented in the
[NNV3.0 README](../README.md) — Option A (`Dockerfile.online`, browser
sign-in) and Option B (root `Dockerfile`, `port@host` licence). Both
flows auto-install AIVL via `mpm`. From inside the resulting container:

```bash
cd code/nnv/examples/NNV3.0/ToolComparison
matlab -nodisplay -batch "run_toolcomparison('mode','smoke')"   # ~1 min
matlab -nodisplay -batch "run_toolcomparison"                    # ~3 h
```

## Install AIVL

Both NNV3.0 Docker flows (Option A `Dockerfile.online`, Option B main
`Dockerfile`) install the AIVL Support Package automatically via `mpm`
(`Deep_Learning_Toolbox_Verification_Library`). Reviewers using either
Docker path do **not** stage anything separately — verify in the smoke
log that the `[AIVL] status: AVAILABLE` line appears.

### Host MATLAB (no Docker)

Users running the experiments directly against a host MATLAB R2025b
install AIVL via the Add-On Explorer:

1. In MATLAB, **Home → Add-Ons → Get Add-Ons**.
2. Search for **"AI Verification Library"** and install.
3. Restart MATLAB. Verify the entry points resolve:

   ```matlab
   which verifyNetworkRobustness        % should resolve to an aivnv dir
   which estimateNetworkOutputBounds
   ```

No further setup is needed — `run_toolcomparison` picks AIVL up
automatically.

### Running without AIVL

If your MATLAB licence doesn't entitle the AI Verification Library, or
you want a faster NNV-only sanity pass, run:

```matlab
run_toolcomparison('mode','default','tools',{'nnv'})
```

The harness works unchanged; AIVL rows are simply absent from
`table_main`. The headline NNV results (paper Tables 5–7 NNV columns)
are unaffected.

## Benchmark catalog

| Benchmark | N inst | Network | Property | NNV (default) | AIVL |
|---|---:|---|---|---|---|
| **acas_xu_p3** | 20 | ACAS Xu FC ReLU (6 layers) | half-space VNNLIB (prop_3) | `approx-star`, `exact-star`, `relax-star-range-50` | `estimate-bounds` |
| **acas_xu_p4** | 20 | ACAS Xu FC ReLU (6 layers) | half-space VNNLIB (prop_4) | `approx-star`, `exact-star`, `relax-star-range-50` | `estimate-bounds` |
| **rl** | 50 | Cartpole + Lunarlander FC (tanh+ReLU) | half-space VNNLIB | `approx-star`, `exact-star`, `relax-star-range-50` | `estimate-bounds` |
| **oval21** | 30 | CIFAR ResNet (3 nets × 10 props) | half-space VNNLIB | `approx-star` | `estimate-bounds` |
| **collins_rul** | 62 | Small 1D Conv (Collins RUL) | half-space VNNLIB | `approx-star` | `estimate-bounds` |
| **mnist_resnet8** | 25 imgs × 4 ε | Custom MNIST ResNet-8 (additionLayer) | argmax robustness (L∞) | `relax-star-area-50` | `deep-poly` (verifyNetworkRobustness) |

In **smoke mode** every benchmark runs one instance × one NNV algorithm
and (when AIVL is on path) one instance × one AIVL algorithm. Smoke
produces 12 result rows (or 6 when AIVL is not installed) in ~1 minute
(the wall is dominated by `parpool` startup + the first AIVL call's
~14 s deep-poly warmup on MNIST-ResNet-8; pure verification work is a
few seconds total). Doubles as an AIVL-setup validation: a
`[AIVL] status: AVAILABLE` line in the terminal/run.log confirms the
Support Package was installed correctly, while a
`[AIVL] status: NOT FOUND` warning surfaces a missing Support Package
(licence not entitled, or `mpm` step skipped) before reviewers waste
time on the full-mode run.

In **default mode** every benchmark runs every algorithm above on its
full instance count. AIVL rows are included.

## Mode reference

| Mode | What runs | Wall-clock | Disk written |
|---|---|---|---|
| `'smoke'` | 1 inst × 1 NNV alg + 1 inst × 1 AIVL alg per benchmark when AIVL is on path (12 rows); 1 inst × 1 NNV alg only otherwise (6 rows) | ~1 min (parpool-startup-dominated) | <1 MB results |
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
- No cgan / TransposedConv2D — AIVL R2025b doesn't support
  `TransposedConv2D` or custom `FunctionLayer`, blocking the cgan
  comparison entirely.

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
│   ├── parse_argmax_vnnlib.m
│   ├── load_mw_network.m
│   └── addpath_shared.m
├── tables/
│   ├── make_table_main.m              consolidated table renderer
│   └── out/                           generated table_main.{tex,txt}
├── results/                           created at runtime
├── logs/                              created by run_all.sh
└── ISSUES.md                          auto-appended summary
```
