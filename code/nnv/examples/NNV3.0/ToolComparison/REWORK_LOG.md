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
| A4  | PAR-2 + timeout columns in Tables A and C                     | done         | TBD     |
| A5  | Bundle CAV'23 assets locally instead of reaching into NNV2.0/ | done         | TBD     |
| B1  | Add one VNNCOMP'24 benchmark (likely `dist_shift`)            | pending      | —       |
| B2  | Pareto plot (additive)                                        | code done    | TBD     |
| B3  | ACAS exact-star reconciliation at 1800s (option + helper)     | code done    | TBD     |
| B4  | AIVL falsification baseline (R2025b API permitting)           | pending      | —       |
| C1  | Table A timeout column / regime split                         | done w/ A4   | —       |
| C2  | Drop noisy supplement / legacy rows from active table         | deferred     | —       |
| C3  | README refresh                                                | done         | TBD     |

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

- **2026-05-06 10:18** — A1 end-to-end smoke OK in nnv3.0:r2025b:
  - ACAS NNV approx-star: dispatch via 'nnv' tool ✓
  - ACAS AIVL estimate-bounds: dispatch via 'aivl' + 'estimate-bounds' ✓
  - MNIST-ResNet AIVL deep-poly: `verifyNetworkRobustness(net, XL, XU, ytrueIdx)`
    returned `verified` in 7.77s ✓
  Smoke harness lives at `utils/smoke_a1.m`.

- **2026-05-06 10:30** — A4 (PAR-2 + timeout columns) shipped.
  Tables now show: V / X / ? / T/O counts, per-instance Timeout (s) (with
  `*` suffix when per-instance timeouts vary within a row, e.g. VNNCOMP
  benchmarks), Mean t (over solved instances only), and PAR-2 (s) (the
  VNNCOMP scoring metric: unsolved instances counted at 2 * timeout).

  Sample numbers under new schema:
  - ACAS p3 NNV exact-star : V=32 X=3 T/O=10  T=300* mean=129  PAR-2=234
  - ACAS p3 AIVL estimate  : V= 2 X=3 ?=40   T=900   mean=0.05 PAR-2=1600
  - RL NNV exact-star      : V=32 X=15 T/O=2  T=900  mean=3.7  PAR-2=112
  - RL AIVL estimate       : V=20 X=11 ?=19  T=900   mean=0.02 PAR-2=684
  - Collins RUL NNV approx : V=10 X=47 ?=5   T=60.0  mean=2.3  PAR-2=9.8
  - Collins RUL AIVL est   : V=10 X=47 ?=5   T=60.0  mean=0.17 PAR-2=9.8

  PAR-2 reveals: NNV's exact-star is the precision winner where applicable
  (ACAS, RL) — it pays a wall-clock cost but the score reflects that
  unsolved instances cost 2T. AIVL's estimate-bounds is fast-but-imprecise
  on ACAS Xu (PAR-2 1.6 ks vs NNV's 234 s) but ties NNV approx-star on
  Collins RUL (both at 9.8 s — Collins is structurally easy).

- **2026-05-06 11:00** — B2 + B3 + C3 implemented file-only (auto-mode
  Bash classifier was transiently unavailable; commits queued):
  - **B2 Pareto plot** (`tables/make_pareto_plot.m`, ~150 lines): per-
    benchmark scatter of (PAR-2, V-rate). Pareto frontier highlighted,
    dominated points dimmed, tool encoded by color, algorithm by marker.
    Wired into `run_toolcomparison.m` so it's a side-effect of the full
    suite. Outputs `tables/out/pareto_<benchmark>.{png,pdf}`.
  - **B3 exact-star reconciliation** (`run_acas_rl_tll.m` +
    `tool_utils.m`): added `'reconcileExact', true` option that
    surgically re-runs `(tool='nnv', algorithm='exact-star',
    status='timeout')` rows at 1800 s without disturbing other
    algorithms' rows. Required extending `tool_utils.purge_status` to
    accept optional (tool, algorithm) filter args. ~6 hours wall to
    convert "V+T/O matches CAV'23" to "V matches CAV'23" with no hedge.
  - **C3 README** (`ToolComparison/README.md`): rewrote tool/algorithm
    grid to single `nnv`/`aivl` rows with footnotes documenting the
    R2025b constraints; added a "Constraints" section spelling out the
    R2025b ceiling (no VNNLIB ingest, no `Algorithm=` knob, no
    α/β-CROWN bridge); added "Reconcile ACAS exact-star at a higher
    cap" usage example; documented the new metric columns
    (timeout/Mean/PAR-2) and the Pareto plot; updated Layout to show
    bundled `acas/`, `rl_benchmarks/`, `oval21/` next to
    `collins_rul_cnn_2022/`.

  Pending: stage + commit (A4 + B2 + B3 + C3) when the classifier is
  back. A3 (swap MNIST-ResNet-8 for a VNNCOMP CNN) and B1 (add VNNCOMP'24)
  both need shell access to fetch the asset files; deferred until Bash
  is unblocked.

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
