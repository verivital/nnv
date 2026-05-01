# Expected Results — NNV3 Tool Paper Reproduction

This file documents the headline numbers each example produces under
its **paper-faithful default configuration**, alongside the verdicts
and timings observed when this artifact was last validated end-to-end.

Use it to eyeball-compare your reproduction:

- **Verdicts** (verified counts, SAT/UNSAT statuses, fair/unfair %)
  reproduce the paper exactly.
- **Timings** depend on hardware (CPU model, GPU presence, MATLAB
  threading); reference numbers below are from a Windows 11 host with
  an NVIDIA RTX 5090 (32 GB, Blackwell), driver 581.95, CUDA 13,
  running the `nnv3.0` Docker image (MATLAB R2024b).

To run the full paper reproduction (including the long ToolComparison
grid; default `bash run_all.sh` uses the ~12 min ToolComparison smoke mode
so the end-to-end suite stays under ~30 min):

```bash
docker run --rm --gpus all \
    -v "$PWD":/home/matlab/nnv \
    -w /home/matlab/nnv/code/nnv/examples/NNV3.0 \
    -e TOOLCOMPARISON_MODE=full \
    nnv3.0 bash run_all.sh
```

After completion, compare your `<Folder>/results/<ts>/` outputs
against the per-example tables below. (See [`README.md`](README.md)
for prerequisites and one-time host setup.)

---

## ModelStar — Paper Table 3 (Appendix B)

**Configuration**: MNIST MLP (5 hidden layers: 1024, 512, 256, 256,
256), 100 test images, single-layer weight perturbation, ReLU
activations, approx-star reachability on CPU. Sweep is over fc_4,
fc_5, fc_6 with the paper's perturbation grid (L∞ as % of layer's
weight range).

| Layer  | Pert. (% of range) | Paper verified | Observed | Paper time / image | Observed time / image |
|--------|-------------------:|---------------:|---------:|-------------------:|----------------------:|
| fc_6   | 0.5 %             | 100 / 100       | 100 / 100 ✓ | 0.18 s | ~0.05 s |
| fc_6   | 1.0 %             | 83 / 100        | 83 / 100  ✓ | 0.15 s | ~0.04 s |
| fc_6   | 1.5 %             | 32 / 100        | 32 / 100  ✓ | 0.11 s | ~0.03 s |
| fc_6   | 2.0 %             | 8 / 100         | 8 / 100   ✓ | 0.07 s | ~0.02 s |
| fc_5   | 0.1 %             | 100 / 100       | 100 / 100 ✓ | 0.72 s | ~0.22 s |
| fc_5   | 0.2 %             | 92 / 100        | 92 / 100  ✓ | 1.04 s | ~0.27 s |
| fc_5   | 0.3 %             | 38 / 100        | 38 / 100  ✓ | 1.37 s | ~0.32 s |
| fc_5   | 0.4 %             | 14 / 100        | 14 / 100  ✓ | 1.41 s | ~0.34 s |
| fc_4   | 0.1 %             | 100 / 100       | 100 / 100 ✓ | 1.94 s | ~0.41 s |
| fc_4   | 0.2 %             | 69 / 100        | 69 / 100  ✓ | 3.33 s | ~0.67 s |
| fc_4   | 0.3 %             | 11 / 100        | 11 / 100  ✓ | 4.67 s | ~0.89 s |
| fc_4   | 0.4 %             | 0 / 100         | 0 / 100   ✓ | 5.75 s | ~1.12 s |

**12 / 12 cells reproduce paper Table 3 exactly.** Observed wall time
for the full sweep: **~70 minutes**.

Output: `ModelStar/results/MNIST_MLP.mat` (in-place population of the
template by `EXPT.save`).

---

## VideoStar — Paper Table 4 (Appendix C)

**Configuration**: ZoomIn-4f, 4-frame MNIST video classifier, 10 test
samples, 30-min timeout per sample, `relax` algorithm.

| ε     | Paper verified / unknown | Observed verified / unknown | Paper avg time | Observed avg time |
|-------|-------------------------:|----------------------------:|---------------:|------------------:|
| 1/255 | 7 / 3                    | 7 / 3 ✓                      | 85.0 s         | ~21 s             |
| 2/255 | 7 / 3                    | 7 / 3 ✓                      | 84.2 s         | ~23 s             |
| 3/255 | 7 / 3                    | 7 / 3 ✓                      | 80.3 s         | ~22 s             |

Observed total wall time: **~12 minutes**.

Outputs: `VideoStar/results/<ts>/eps=1_255.csv`,
`eps=2_255.csv`, `eps=3_255.csv`.

---

## GNNV — Paper Figure 5 (Appendix D)

**Configuration**: PowerFlow on IEEE 24-bus, 10 graphs
(filtered to the "safe" GT-voltage subset), node-only perturbation,
4-ε sweep.

| Architecture | ε = 1e-5 | ε = 1e-4 | ε = 1e-3 | ε = 1e-2 |
|--------------|---------:|---------:|---------:|---------:|
| GCN          | 6.9 %    | 6.2 %    | 5.4 %    | 0.0 %    |
| SAGE         | 99.2 %   | 99.2 %   | 99.2 %   | 99.2 %   |
| GINE-Conv    | 87.7 %   | 87.7 %   | 86.2 %   | 27.7 %   |

**Paper's qualitative finding reproduces**: SAGE is the most robust
across all perturbations, GCN is the least, GINE-Conv is intermediate.
Observed total wall time: **~3.5 minutes**.

Output: `GNNV/results/gnn_<ts>/gnn_results.csv` (one row per arch × task
× mode × grid × ε).

---

## ProbVer — Paper Table 5 (Appendix E)

**Configuration**: TinyYOLO from VNN-COMP 2023, 3 randomly selected
VNNLIB properties from the 72-instance benchmark, ε = 1/255, cp-star
reachability, surrogate model with 99.9 % coverage / 99.9 % confidence.
Random-instance seed: 42.

| Property (instance #) | Paper verdict | Observed verdict | Paper time (CPU) | Observed time (CPU) |
|-----------------------|---------------|------------------|-----------------:|--------------------:|
| Prop 27               | (varies)      | UNSAT ✓           | 528–602 s        | ~163 s              |
| Prop 52               | (varies)      | UNSAT ✓           | 528–602 s        | ~163 s              |
| Prop 68               | (varies)      | UNSAT ✓           | 528–602 s        | ~163 s              |

The paper randomly selected 3 properties without disclosing the seed,
so the **specific instance indices may differ**. The expected
qualitative pattern matches: **all 3 properties verify as UNSAT**.
Observed total wall time: **~9 minutes** (519 s wall, of which ~489 s
is the three per-instance MATLAB processes — cf.\ `~163 s` per instance
in the table — plus ~30 s of inter-process startup).

The 528–602 s paper figures came from CPU-only runs on a machine with
a slower single-thread profile and full TinyYOLO surrogate training each
time. Observed ~163 s reflects the same code path in this artifact's
container with a faster CPU and warm BLAS / MEX caches; the verdicts are
unchanged. Either bracket is acceptable for repeatability.

Output: `ProbVer/results/<ts>/results_summary.csv`.

---

## FairNNV — Paper Table 6 (Appendix F)

**Configuration**: UCI Adult-Income, 100 test observations, two
classifiers (AC-1 small: 13 → 16 → 8 → 2; AC-3 medium: 13 → 50 → 2),
seed 500.

### Counterfactual fairness (ε = 0)

| Model      | Paper VF | Observed VF | Paper time / sample | Observed time / sample |
|------------|---------:|------------:|--------------------:|-----------------------:|
| AC-1 small | 89 %     | 89 % ✓       | 0.0211 s            | ~0.009 s               |
| AC-3 med   | 87 %     | 87 % ✓       | 0.0190 s            | ~0.004 s               |

### Individual fairness — AC-1 (small, 16-8)

| ε    | Paper VF | Observed VF | Paper time / sample | Observed time / sample |
|------|---------:|------------:|--------------------:|-----------------------:|
| 0.01 | 87 %     | 87 % ✓       | 0.0235 s            | ~0.006 s               |
| 0.02 | 84 %     | 84 % ✓       | 0.0278 s            | ~0.006 s               |
| 0.03 | 81 %     | 81 % ✓       | 0.0374 s            | ~0.008 s               |
| 0.05 | 69 %     | 69 % ✓       | 0.0557 s            | ~0.012 s               |
| 0.07 | 50 %     | 50 % ✓       | 0.0882 s            | ~0.020 s               |
| 0.10 | 22 %     | 22 % ✓       | 0.1546 s            | ~0.033 s               |

### Individual fairness — AC-3 (medium, 50)

| ε    | Paper VF | Observed VF | Paper time / sample | Observed time / sample |
|------|---------:|------------:|--------------------:|-----------------------:|
| 0.01 | 86 %     | 86 % ✓       | 0.0671 s            | ~0.015 s               |
| 0.02 | 84 %     | 84 % ✓       | 0.1531 s            | ~0.032 s               |
| 0.03 | 82 %     | 82 % ✓       | 0.2648 s            | ~0.052 s               |
| 0.05 | 71 %     | 71 % ✓       | 0.5703 s            | ~0.112 s               |
| 0.07 | 50 %     | 50 % ✓       | 1.1121 s            | ~0.211 s               |
| 0.10 | 27 %     | 27 % ✓       | 2.7662 s            | ~0.539 s               |

**All VF percentages reproduce paper Table 6 exactly.** Per-sample
timings on are based on hardware. Observed total wall time: **~2 minutes**.

Outputs (under `FairNNV/results/<ts>/`):

- `counterfactual_<ts>.csv`, `individual_<ts>.csv`, `timing_<ts>.csv`
- `counterfactual_table.tex`, `timing_table.tex` (LaTeX tables matching the paper)
- `individual_fairness_combined.{png,pdf}` (Figure 7)

---

## ToolComparison — Paper Tables 5, 6, and 7

**Configuration**: NNV vs MathWorks AI Verification Library (AIVL) on
five VNNLIB benchmarks (ACAS Xu p3/p4, RL controllers, OVAL21, Collins
RUL CNN) plus MNIST-ResNet-8 argmax robustness. AIVL Support Package
required for the `mw_*` rows; pass `'tools',{'nnv'}` to skip. Two modes:

- `'mode','smoke'` — 5 ACAS networks + 5 MNIST images, NNV-only,
  `approx-star` + `relax-star-area-50`. Used by `run_all.sh` to keep
  the end-to-end suite under ~30 min.
- `'mode','full'` (default) — full per-benchmark NNV grid + AIVL.
  Renders the paper's Tables 5, 6, and 7.

### Paper Table 5 (FC VNNLIB benchmarks)

`tables/out/table_A.{tex,txt}`, full-mode only. Headline NNV cell per
benchmark vs AIVL `estimate-bounds`:

| Benchmark           | NNV best                | AIVL                | Headline finding |
|---------------------|-------------------------|---------------------|------------------|
| acas_p3 (45)        | exact-star V=32, X=3, T=10 | V=2, X=3, ?=40   | NNV verifies 16× more |
| acas_p4 (45)        | exact-star V=39, X=3, T=3  | V=0, X=1, ?=44   | NNV verifies 39 of 45 |
| rl   (50)           | approx-star V=32, X=14     | V=20, X=11, ?=19 | NNV verifies 60% more |

`tables/out/sanity_report.txt` confirms ACAS p3/p4 match CAV'23
exact-star (`Verified + Timeout = 42` exactly).

### Paper Table 6 (CNN VNNLIB benchmarks)

Same `tables/out/table_A.{tex,txt}` (rendered together with Table 5):

| Benchmark           | NNV best                 | AIVL              | Headline finding |
|---------------------|--------------------------|-------------------|------------------|
| oval21 (30)         | approx-star X=10 (5 timeouts) | X=9          | NNV +1 counterexample |
| collins_rul (62)    | approx-star V=10, X=47, ?=5 | V=10, X=47, ?=5 | identical verdicts |

### Paper Table 7 (MNIST-ResNet-8)

`tables/out/table_C.{tex,txt}`, full-mode only. Both tools verify all 50
images per ε; NNV's `relax-star-area-50` is faster:

| ε      | AIVL deep-poly mean t | NNV mean t | NNV speedup |
|--------|----------------------:|-----------:|------------:|
| 1/255  | 45.7 s                | 6.2 s      | 7.4× |
| 2/255  | 43.0 s                | 8.3 s      | 5.2× |
| 4/255  | 42.9 s                | 13.7 s     | 3.1× |
| 8/255  | 41.4 s                | 38.6 s     | 1.1× |

### Smoke mode (used by `run_all.sh`)

Smoke mode produces a much smaller subset (5 ACAS networks × 2
algorithms + 5 MNIST images × 2 algorithms, NNV-only) — sanity check,
not paper-grade. Expected wall time **~12 minutes**.

**Full-mode total wall time**: ~3–5 h on a 4-core CPU (driven by
`exact-star` on ACAS Xu, ≤900 s/instance). Set
`TOOLCOMPARISON_MODE=full bash run_all.sh` to opt in.

Outputs (under `ToolComparison/`):

- `acas_rl_tll/results/results_<benchmark>.{mat,csv}`
- `mnist_resnet/results/expC_<model>.{mat,csv}`
- `tables/out/{table_A,table_C}.{tex,txt}` (paper Tables 5/6, and Table 7)
- `tables/out/sanity_report.txt`

---

## Determinism notes

- Each example seeds its RNG (FairNNV = 500, ProbVer = 42; GNNV /
  ModelStar / VideoStar use deterministic test-set indices).
  Verdicts are exact-reproducible.
- Timings vary with hardware. Don't expect the *exact* numbers in the
  observed columns above — use them as an order-of-magnitude reference.

## Total wall-clock for the full sweep

On the reference host:

| Example | Wall time |
|---------|----------:|
| FairNNV | ~2 min    |
| GNNV    | ~3.5 min  |
| ProbVer | ~9 min    |
| VideoStar | ~12 min |
| ModelStar | ~70 min |
| ToolComparison (smoke) | ~12 min |
| ToolComparison (full)  | ~3–5 h (opt-in via `TOOLCOMPARISON_MODE=full`) |
| **Total (smoke ToolComparison)** | **~107 min (~1.8 h)** |
| **Total (full ToolComparison)**  | **~5–7 h** |
