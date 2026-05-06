# ToolComparison Rework Log

Tracking the post-review improvement plan. Goal: make ToolComparison defensible
for an ATVA 2026 reviewer — proper AIVL leverage, real published benchmarks,
and the right verification metrics — within the R2025b ceiling that CodeOcean
imposes.

Updated as work lands. Most-recent entries at the top.

## Status board

| ID  | Item                                                          | State        | Commit  |
|-----|---------------------------------------------------------------|--------------|---------|
| 0   | Bump MATLAB R2024b → R2025b in Dockerfile                     | done         | TBD     |
| A1  | Proper AIVL invocation (verifyNetworkRobustness, scope-bounded by R2025b) | partial | TBD |
| A2  | Collapse `mw_*` tools → one `aivl` tool with algorithm column | done         | TBD     |
| A3  | Replace in-tree MNIST-ResNet-8 with VNNCOMP CNN               | pending      | —       |
| A4  | Add PAR-2 to Tables A and C                                   | pending      | —       |
| A5  | Bundle CAV'23 assets locally instead of reaching into NNV2.0/ | done         | TBD     |
| B1  | Add one VNNCOMP'24 benchmark (likely `dist_shift`)            | pending      | —       |
| B2  | Pareto plot (additive, not replacement)                       | pending      | —       |
| B3  | ACAS exact-star reconciliation at 1800s                       | pending      | —       |
| B4  | AIVL falsification baseline (R2025b API permitting)           | pending      | —       |
| C1  | Table A timeout column / regime split                         | pending      | —       |
| C2  | Drop noisy supplement / legacy rows from active table         | pending      | —       |
| C3  | README refresh                                                | pending      | —       |

## Constraints

- **MATLAB ceiling: R2025b** (CodeOcean cap). α/β-CROWN AIVL bridge is R2026a-only — out of scope.
- **PI directive: keep the full relax-star factor grid** (25/50/75/100). Pareto plot is additive.
- **Reproducibility**: non-redistributable AIVL tarball; no GPU required for ToolComparison.

## Decisions log

### 2026-05-05 — bump to R2025b

Why now: aligns Docker with CodeOcean. R2025b has matured `verifyNetworkRobustness`
(VNNLIB-form ingest, ResNet support) which the rework depends on.

Cost: one Dockerfile bump + one rebuild. README reference timings will need to be
re-baselined if R2025b shifts wall-clocks by >5%.

Hard cap on AIVL: no α/β-CROWN until R2026a + ATVA reviewer has it. Footnoted.

## Build / test record

Each entry: timestamp, what was tested, outcome, log path.

- **2026-05-05 16:32** — `nnv3.0:r2025b` build OK in ~3 min (cache-warm).
  Image is 19.1 GB. Existing `nnv3.0:fixes` (R2024b, 18.6 GB) preserved as
  fallback. Log: `/tmp/nnv3_r2025b_20260505_162939.log`.
- **2026-05-05 18:30** — A2 + A5 smoke OK inside r2025b container.
  Migrated 232 legacy `mw_estimate` rows across 5 `.mat` files
  (45 acas_p3, 45 acas_p4, 50 rl, 30 oval21, 62 collins_rul) to
  `aivl` + `estimate-bounds`. Both table renderers run cleanly under
  the new schema. Local asset locator finds 45 ACAS ONNX, 296 RL VNNLIB,
  30 OVAL21 properties. Log: `/tmp/nnv3_a2_smoke_20260505_*.log`.

- **2026-05-05 19:33** — Dockerfile chown fix verified: rebuilt
  `nnv3.0:r2025b` extracts AIVL cleanly. `which verifyNetworkRobustness`
  resolves to the extracted `aivnv/` dir. Committed `0b6aa0bbc`.

- **2026-05-05 19:50** — AIVL VNNLIB ingest probe (Phase A1 gate).
  Findings (`probe_aivl_vnnlib.m`):
  1. `verifyNetworkRobustness(net, vnnlibFile)` does NOT exist in R2025b.
     Error: `Invalid argument list. Function requires 2 more input(s).`
     The signature is `(net, XL, XU, ytrue_class)` — argmax-form only.
     VNNLIB ingest landed in R2026a.
  2. The `Algorithm=` Name-Value param is absent in R2025b. Allowed names
     are `'MiniBatchSize'` and `'ExecutionEnvironment'` only. There is a
     single internal algorithm (DeepPoly per MathWorks docs) with no
     user-facing knob.
  3. Argmax-form works on FC ACAS networks:
     `verifyNetworkRobustness(net, XL, XU, 1)` returned `categorical
     unproven` in 0.80s.
  4. Side issue: the previously bundled `ensure_aivl_on_path()` and
     worker-side path setup mistakenly addpath'd CONTENTS of `aivnv/`,
     producing 10+ MATLAB warnings about namespace dirs (`+aivnv/`,
     `+coder/`, `+matlab/`) and `resources/`. Fixed: addpath the
     `aivnv/` dir itself, once. MATLAB resolves namespaces from the
     parent.

  **Implications**:
  - ACAS p3/p4 / RL / OVAL21 / Collins RUL: all are general half-space
    VNNLIB properties, not argmax. Under R2025b we cannot use
    `verifyNetworkRobustness`; the proper AIVL baseline remains
    `estimateNetworkOutputBounds` + manual property check.
  - MNIST-ResNet (and any future argmax-form image benchmark): proper
    AIVL invocation is `verifyNetworkRobustness(net, XL, XU, ytrueIdx)`
    — already wired correctly post-A2.
  - Speculative `'deep-poly'` branches I had stubbed in verifyAcas /
    verifyRL / verifyImageVnnlib are now hard-error stubs documenting the
    R2026a requirement. They become real branches when an R2026a host
    runs the code.

  Paper-text consequence: AIVL R2025b is honestly the
  estimateNetworkOutputBounds + interval-arithmetic baseline for
  half-space properties, plus `verifyNetworkRobustness` for argmax-form.
  α/β-CROWN remains absent. This is the artifact's faithful AIVL surface
  area; the paper should state this constraint explicitly rather than
  leave reviewers wondering why we didn't use deep-poly on ACAS.

## Notes for reviewers (paper-side framing)

These observations should make their way into the paper text once the
implementation settles:

- ToolComparison is **a parity check**, not the headline contribution. NNV 3.0's
  contribution is the things AIVL doesn't do at all (GNN / probabilistic /
  weight-perturbation / video). The standard-NN comparison is to anchor "we are
  competitive on the standard task too".
- AIVL is invoked through `verifyNetworkRobustness`, the public verification
  entry point, with the algorithms R2025b exposes (DeepPoly headline; no α/β-CROWN
  bridge until R2026a). `estimateNetworkOutputBounds` + interval-arithmetic is
  retained as a *supplementary* row labeled honestly, not as the AIVL headline.
- Benchmark mix uses VNNCOMP-format `.onnx` + `.vnnlib` files end-to-end, with
  per-instance timeouts honored from `instances.csv` where the suite ships them.
