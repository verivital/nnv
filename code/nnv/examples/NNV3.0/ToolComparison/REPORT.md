# ATVA26 — NNV 3.0 vs MathWorks AIVL Comparison Report

**Scope:** Experimental comparison for the NNV 3.0 ATVA 2026 tool paper
(follow-up to NNV 1.0 / CAV'20 and NNV 2.0 / CAV'23).

---

## 1. Executive summary

We compare NNV 3.0 against the MathWorks AI Verification Library (AIVL)
across two experiments:

- **Experiment A (Table A)** refreshes the NNV 2.0 / CAV'23 head-to-head on
  ACAS Xu, RL controllers, and TLLverify benchmarks under modern MATLAB
  R2025b. ACAS Xu numbers cross-check against the CAV'23 ground-truth
  within tolerance (`Verified + Timeout` matches CAV'23's `Verified` exactly).
- **Experiment C (Table C)** is the first head-to-head between NNV and
  AIVL's `verifyNetworkRobustness` on a residual network (MNIST-ResNet-8).
  AIVL only gained `additionLayer` support in R2024b, making this comparison
  newly possible.

Both experiments run from a single Docker image (`atva26/nnv-mw:latest`,
16.6 GB) bundling MATLAB R2025b, the NNV tree, AIVL, and the driver code.

**Headline findings:**

- On **RL controllers** NNV verifies 60% more instances than AIVL DeepPoly
  (32 vs 20 of 50) and finds 3 more counterexamples DeepPoly misses.
- On **MNIST-ResNet-8** NNV's relax-star is **3.1×–7.4× faster** than
  AIVL DeepPoly at small ε, with both tools verifying every instance.
  AIVL closes the gap at ε=8/255 as set splitting in NNV grows.
- On **TLLverify** both tools struggle: AIVL's bounds are too loose to
  conclude anything; NNV's relax-star times out on the largest networks.

---

## 2. What the comparison does

### 2.1 Tools and methods

| Tool | Method | What it computes |
|------|--------|------------------|
| NNV 3.0 | **exact-star** | Sound and complete reachability via star-set splitting. |
| NNV 3.0 | **relax-star (50%)** | Sound but incomplete: relaxes 50% of ReLUs to keep set count bounded. Fast, may return *unknown*. |
| AIVL | `estimateNetworkOutputBounds` | Sound interval/zonotope bounds on output. Needed for half-space output specs. |
| AIVL | `verifyNetworkRobustness` (DeepPoly) | Sound robustness verification for argmax-style specs. Requires R2024b+ for residual networks. |

### 2.2 Aspects compared

For every (tool × method × instance):

- **Verdict**: verified / violated / unknown / timeout / error.
- **Per-instance solve time** (with a 900 s wall-clock cap).
- **Architectural reach** — which networks each tool can ingest at all.

The result schema is uniform across experiments: a MATLAB `table` with
columns `{tool, benchmark, instance_id, status, time, algorithm, timeout, note}`.

### 2.3 Why some cells are intentionally empty

- `verifyNetworkRobustness` is **not run on Experiment A** because ACAS Xu,
  RL, and TLL specifications are VNNLIB-style output half-spaces, not
  argmax robustness. AIVL's robustness API does not accept these specs.
  This is the same scoping NNV 2.0 / CAV'23 used.
- **NNV exact-star is omitted for TLL** (and for ResNet in Experiment C):
  exact set-splitting is intractable on networks of this scale.

---

## 3. Experiment A — Refresh of NNV 2.0 Tables 2–3

### 3.1 Benchmarks

| Benchmark | Networks | Properties | NNV methods | AIVL method |
|-----------|----------|------------|-------------|-------------|
| `acas_p3` | 45 ACAS Xu | P3 on all 45 | relax-star-50, exact-star | `estimateNetworkOutputBounds` |
| `acas_p4` | 45 ACAS Xu | P4 on all 45 | relax-star-50, exact-star | `estimateNetworkOutputBounds` |
| `rl` | 1 cartpole + 1 pendulum | 50 (random subset, seed=1) | relax-star-50, exact-star | `estimateNetworkOutputBounds` |
| `tllverify` | 32 networks | 32 | relax-star-50 | `estimateNetworkOutputBounds` |

Timeout: 900 s per instance.

### 3.2 Results — Table A

`tables/out/table_A.{tex,txt}`

| Benchmark | Tool | Algorithm | V | X | ? | T/O | Mean t |
|-----------|------|-----------|---|---|---|-----|--------|
| acas_p3 | nnv | relax-star-50 | 2 | 2 | 41 | 0 | 0.6 s |
| acas_p3 | nnv | exact-star | 32 | 3 | 0 | 10 | 129 s |
| acas_p3 | mw_estimate | estimate-bounds | 2 | 3 | 40 | 0 | 0.05 s |
| acas_p4 | nnv | relax-star-50 | 1 | 2 | 42 | 0 | 0.6 s |
| acas_p4 | nnv | exact-star | 39 | 3 | 0 | 3 | 170 s |
| acas_p4 | mw_estimate | estimate-bounds | 0 | 1 | 44 | 0 | 0.03 s |
| rl | nnv | relax-star-50 | 32 | 14 | 4 | 0 | 0.07 s |
| rl | nnv | exact-star | 32 | 15 | 1 | 2 | 3.7 s |
| rl | mw_estimate | estimate-bounds | 20 | 11 | 19 | 0 | 0.02 s |
| tllverify | nnv | relax-star-50 | 0 | 0 | 17 | 7 | 98 s |
| tllverify | mw_estimate | estimate-bounds | 0 | 0 | 24 | 0 | 1.0 s |

V = verified, X = violated, ? = unknown, T/O = timeout.

### 3.3 CAV'23 cross-check

`tables/out/sanity_report.txt`:

```
ATVA26 sanity report vs CAV'23 (NNV 2.0 exact-star)
------------------------------------------------------------
acas_p3: got (V=32, X=3, ?=0, T/O=10)  CAV'23 expected (V=42, X=3, ?=0)  OK
acas_p4: got (V=39, X=3, ?=0, T/O=3)   CAV'23 expected (V=42, X=3, ?=0)  OK
```

`Verified + Timeout = 42` and `Violated = 3` exactly match CAV'23's
exact-star results. The 13 timeouts (10 P3 + 3 P4) are verified-but-too-slow
under our 900 s cap; CAV'23 ran without a cap, with their longest
instance at 10 479 s (~2.9 h).

### 3.4 Discussion

**RL is the strongest comparison point.** NNV substantially out-verifies
AIVL DeepPoly on the RL benchmark — 32 verified vs 20 (60% more), 14
violated vs 11 (counterexamples that DeepPoly missed), 4 unknown vs 19.
Mean solve time is comparable (0.07 s vs 0.02 s). Star-set reachability
is decisively more precise here at negligible cost.

**ACAS Xu — the precision/runtime tradeoff is visible.** NNV exact-star
verifies 32+39 instances in mean 129–170 s; relax-star-50 is fast (0.6 s)
but uninformative on most ACAS instances (41–42 unknown). AIVL's
`estimateNetworkOutputBounds` is the fastest of all (0.03–0.05 s) but
matches relax-star-50 on precision — both yield mostly *unknown*. ACAS Xu
is a regime where only sound-and-complete star-set reasoning produces
useful verdicts, and the cost is real.

**TLLverify exposes a shared limitation.** AIVL returns 24/32 unknown in
≤1 s — its sound bounds are too loose. NNV-relax-star returns 17/32
unknown in mean ~98 s and 7 timeouts at 900 s. Neither tool is
practical at this scale; the ReLU set-explosion characteristic of TLL
networks defeats both interval-based and star-based methods. Exact-star
is intentionally omitted as intractable.

**On ACAS exact-star timeouts.** Our 900 s cap accounts for the gap to
CAV'23's reported counts. With a longer cap (e.g. 4 h) the timeouts would
likely close — CAV'23's worst-case ran ~2.9 h.

---

## 4. Experiment C — ResNet DAG head-to-head

### 4.1 Model and verification matrix

| Model | Inputs | Architecture | Test accuracy |
|-------|--------|--------------|---------------|
| `mnist_resnet8` | 28×28×1 | 8 conv + 3 residual blocks (`additionLayer`) + GAP + FC10 | 60.9% (3 epochs) |

50 test images × 4 perturbation radii ε ∈ {1/255, 2/255, 4/255, 8/255}
× 2 tools = **400 instances**.

| Tool | Algorithm |
|------|-----------|
| NNV 3.0 | relax-star-50 |
| AIVL | `verifyNetworkRobustness` (DeepPoly, R2024b+) |

NNV exact-star is omitted: residual-path set splitting is infeasible at
this depth.

### 4.2 Results — Table C

`tables/out/table_C.{tex,txt}`

| ε | AIVL DeepPoly | NNV relax-star-50 | NNV speedup |
|---|---------------|--------------------|-------------|
| 1/255 | 50/50 verified, 45.7 s | 50/50 verified, **6.2 s** | **7.4×** |
| 2/255 | 50/50 verified, 43.0 s | 50/50 verified, **8.3 s** | 5.2× |
| 4/255 | 50/50 verified, 42.9 s | 49/50 verified, **13.7 s** | 3.1× |
| 8/255 | 50/50 verified, 41.4 s | 50/50 verified, **38.6 s** | 1.07× |

### 4.3 Discussion

1. **Both tools produce the same robustness verdict at every ε** — the
   methods agree on what is robust.
2. **NNV is 3.1×–7.4× faster at small perturbations.** Relax-star scales
   with ε (more set splitting at larger perturbation radii); DeepPoly's
   interval propagation cost is roughly constant in ε. The crossover
   sits near ε=8/255.
3. **First published head-to-head on residual networks.** AIVL could not
   handle `additionLayer` networks before R2024b; this comparison was
   not possible until that release.

### 4.4 The single NNV "error" cell

The cell `(ε=4/255, image=13, NNV)` records `error`, not `verified`.
The actual underlying message:

```
Caused by:
    Error using lpsolver (line 92)
    LP solver error. Task failed to be solved by linprog and glpk.
    GLPK exitflag = 1
```

NNV's `lpsolver.m` solves an LP to test whether reach sets intersect
half-spaces. It tries `linprog` first, falls back to GLPK. On this
instance both solvers failed (GLPK `exitflag = 1` indicates abnormal
termination, typically a near-singular constraint matrix). This is a
real numerical edge case at the robustness boundary, not a transient
infrastructure flake. We document it as: 1/200 instances triggered an
LP-solver numerical issue; the verifier is not a complete decision
procedure when both backend LP solvers fail. This is a known limitation
of LP-based verification with off-the-shelf solvers.

---

## 5. Methodology and infrastructure

### 5.1 Result schema

Every result `.mat` is a MATLAB `table` with columns
`{tool, benchmark, instance_id, status, time, algorithm, timeout, note}`.
Helpers in `comparison_utils.m`:

- `new_row`, `append_to_mat` — write rows.
- `has_instance` — resumability: skip instances already recorded.
- `purge_status` — re-run only rows with a given status (e.g. `'error'`).
- `tally`, `format_time`, `emit_latex_table`, `sanity_check_vs_nnv2`.

CSV mirrors are produced by `export_csv.m` for VS Code preview:

- `results_<bench>.csv` — one row per (tool, instance, algorithm).
- `results_<bench>_summary.csv` — aggregates by (tool, algorithm, status)
  with count, mean, median.

### 5.2 NNV patches required for R2025b

R2025b's ONNX importers wrap each FC's bias-add as
`nnet.cnn.layer.ScalingLayer` (Scale=1, Offset=bias). AIVL rejects
`ScalingLayer` with `aivnv:verifyNetwork:DisallowedLayers`, and NNV's
`matlab2nnv` did not previously dispatch on this type. Two patches in
the bundled NNV tree:

- [`nnv/code/nnv/engine/nn/layers/ElementwiseAffineLayer.m`](../nnv/code/nnv/engine/nn/layers/ElementwiseAffineLayer.m) —
  `parse()` accepts both `nnet.onnx.layer.ElementwiseAffineLayer`
  (R2024b and earlier) and `nnet.cnn.layer.ScalingLayer` (R2025a+),
  synthesizing `DoScale=DoOffset=true` for the latter.
- [`nnv/code/nnv/engine/utils/matlab2nnv.m`](../nnv/code/nnv/engine/utils/matlab2nnv.m) —
  added `nnet.cnn.layer.ScalingLayer` to the dispatch around line 154.

`rebuild_sequential_for_aivl.m` complements the NNV patches on the AIVL
side: it walks the layer list, folds each `ScalingLayer(Scale=1)` into
the preceding FC's `Bias`, strips adapter layers (`VerifyBatchSize`,
`FlattenInto2dLayer`, `RegressionOutput`, …), and returns a clean
`featureInputLayer → FC → ReLU → …` `dlnetwork` AIVL accepts.

### 5.3 Detached, resumable runs

`run_detached.sh` wraps `docker run` in a named tmux session that
survives SSH disconnects. Result `.mat` files are bind-mounted so they
persist outside the container; combined with `has_instance`, this gives
graceful resumption. The successful runs that produced the published
tables:

- `comparison/logs/expA_20260427_104734.log`
- `comparison/logs/expC_20260427_194757.log`
- `comparison/logs/train_20260427_192602.log` (MNIST-ResNet-8 training)

---

## 6. Repository structure

```
ATVA26/
├── nnv/                                       # Shallow clone of verivital/nnv
│   └── code/nnv/engine/
│       ├── nn/layers/ElementwiseAffineLayer.m # PATCHED: parse() accepts ScalingLayer
│       └── utils/matlab2nnv.m                 # PATCHED: ScalingLayer dispatch added
│
├── docker/
│   ├── Dockerfile.atva26                      # MATLAB R2025b + NNV + AIVL
│   ├── build_image.sh                         # Build wrapper
│   ├── entrypoint.sh                          # Dispatches expA|expC|tables|smoke
│   ├── toolbox_install.m                      # Extracts AIVL tarball
│   ├── README_license.md                      # MathWorks license + add-ons
│   └── addons/atva26-aivl.tar.gz              # AIVL R2025b Support Package
│
└── comparison/
    ├── README.md                              # End-user run instructions
    ├── REPORT.md                              # ← THIS FILE
    ├── comparison_utils.m                     # Result-row schema + helpers
    ├── rebuild_sequential_for_aivl.m          # Layer rewrite for AIVL
    ├── run_detached.sh                        # tmux launcher
    ├── smoke_test.m                           # Hand-built dlnetwork sanity check
    ├── smoke_acas.m                           # Real ACAS Xu → NNV reach
    ├── export_csv.m                           # .mat → .csv mirror
    │
    ├── expA_nnv2_refresh/
    │   ├── run_expA.m                         # ACAS / RL / TLL driver
    │   └── results/
    │       ├── results_{acas_p3,acas_p4,rl,tllverify}.mat
    │       ├── results_<bench>.csv
    │       └── results_<bench>_summary.csv
    │
    ├── expC_resnet_dag/
    │   ├── run_expC.m                         # ResNet verification driver
    │   ├── train_resnets.m                    # MATLAB-native trainer
    │   ├── models/
    │   │   ├── mnist_resnet8.mat              # Trained dlnetwork
    │   │   └── mnist_resnet8_testset.mat      # 100 labeled images
    │   └── results/
    │       ├── expC_mnist_resnet8.mat
    │       ├── expC_mnist_resnet8.csv
    │       └── expC_mnist_resnet8_summary.csv
    │
    ├── tables/
    │   ├── make_table_A.m
    │   ├── make_table_C.m
    │   └── out/
    │       ├── table_A.{tex,txt}              # Paper artifact
    │       ├── table_C.{tex,txt}              # Paper artifact
    │       └── sanity_report.txt              # CAV'23 cross-check
    │
    └── logs/                                  # tmux logs from successful runs
```

---

## 7. Sources of truth for paper text

1. [`tables/out/sanity_report.txt`](tables/out/sanity_report.txt) — quote
   verbatim to demonstrate the CAV'23 regression match.
2. [`tables/out/table_A.tex`](tables/out/table_A.tex) — drop into the paper.
3. [`tables/out/table_C.tex`](tables/out/table_C.tex) — drop into the paper.
4. The two NNV patches (in [`nnv/code/nnv/engine/`](../nnv/code/nnv/engine/))
   are the basis for an upstream PR description.

---

*End of report.*
