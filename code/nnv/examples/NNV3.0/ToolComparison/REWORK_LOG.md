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
| A1  | `mw_estimate` → `verifyNetworkRobustness(net, vnnlib, ...)`   | pending      | —       |
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
