# NNV — Post-Merge Roadmap (after PR #290 lands)

*2026-06-10 (Opus 4.8). Forward plan for what to build once the `transformer` PR #290 is merged. Grounds
on [`PR_MERGE_READINESS_PLAN.md`](PR_MERGE_READINESS_PLAN.md) §4 (F1–F5), the VNN-COMP status (§2), and
[`CR_TRANSFORMER_CLAUDE_FINAL.md`](CR_TRANSFORMER_CLAUDE_FINAL.md) residuals, expanded to cover the
adjacent workstreams (SLM, cloud compute, CI modernization, pre-existing subsystem bugs).*

## Where we are at merge

PR #290 ships **sound single-token / FC-simulated transformer verification + fail-loud multi-token guards
+ ONNX/MATLAB import-correctness fixes + a real CI gate + a tbxmanager toolbox cache**. Every reach is a
sound over-approximation or a fail-loud refusal; **no path certifies not-robust/unsafe from an
over-approximation**. The headline gap that remains is **trustworthy multi-token / real-attention ViT
verification** — today it is sound *by refusal*, not by capability. Soundness is **empirical
(MC-containment), not formal proof**.

---

## Tier 0 — Immediately after merge (housekeeping, low effort, time-sensitive)

- **T0.1 — Branch protection on the matrix CI.** Make REQUIRED: `test_transformer_soundness_regression`,
  `test_vnncomp25_regression`, and the "missing/crashed shard reds the build" gate (`ci_report.py`).
  **Documented rule:** no soundness test may ever enter `ci_allowed_failures.txt` (that would defeat the
  gate; the allow-list is only for env-flaky tutorials + explicitly-WIP demos).
- **T0.2 — CI Node 24 migration (HARD DEADLINE).** GitHub forces Node 24 as the default on **2026-06-16**
  and removes Node 20 on **2026-09-16**. The current actions (`setup-matlab@v2` on the legacy workflows,
  plus `cache@v4`/`checkout@v4`/`upload-artifact@v4`/`run-command@v2`) emit the Node-20 deprecation
  warning. Bump legacy `ci.yml` + `regression-tests.yml` `setup-matlab@v2 → @v3` (also unifies the toolbox
  set with the matrix), and confirm the rest are Node-24-ready. Cleanest end-state: **retire the legacy
  single-job `ci.yml`/`regression-tests.yml` in favor of the matrix-sharded CI** (now toolbox-cache-
  accelerated and flake-resistant), keeping only what the matrix doesn't cover.
- **T0.3 — Resolve the [42] exact-star policy.** Current behavior downgrades not-robust→unknown for
  over-approximate layers under `exact-star`. If the maintainer prefers **hard-rejecting** `exact-star` on
  non-exact layers instead, that's a one-line policy change in `NN.m`'s gate. Pick one and document it.
- **T0.4 — Self-contain the Python ONNX importer.** Copy `tools/onnx2nnv_python/{onnx2nnv.py,
  xvalidate.py}` into `code/nnv/tools/onnx2nnv_python/` (the `load_nnv_from_mat.m` header already points
  there). Keeps the VNN-COMP import path reproducible inside the repo.
- **T0.5 — Model-blob hygiene.** Prefer a `download_models.m` + checksum (the `SST2/download_sst2.m`
  pattern) over committing the multi-MB `mnist_vit_*` / `sentiment_*` `.mat`/`.onnx` blobs upstream. If any
  are kept in-tree, call them out explicitly in the PR/release description.

---

## Tier 1 — Headline capability: sound transformer verification

- **T1.1 / F1 (HIGH) — Sound multi-token attention bound.** Replace the `multiTokenUnsound` refusal with a
  real bound. Cheapest sound version (valid because attention weights live in the probability simplex):
  per output coordinate `out_(i,d) ∈ [ min_j V_lb(j,d), max_j V_ub(j,d) ]` (value-hull bounding box,
  ignores Q/K). Requires the reach API to **preserve `[seq_len, dim]`** instead of flattening. **Verify by
  MC for seq_len ∈ {2,4,8}** (0 containment violations). Refs: Shi 2020 (arXiv:2002.06622), Wei 2023
  (arXiv:2303.01713), DeepT/Bonaert 2021, Vertex-Softmax 2026 (arXiv:2605.10974).
- **T1.2 / F2 (HIGH) — Real-attention ViT pipeline.** With F1 + the EmbedDim fix already in the PR, make
  `verify_mnist_vit_attention.m` produce **sound** verdicts, then **un-mark it EXPERIMENTAL**; revisit the
  `single`-precision conversion warning. This is the demo that closes the loop on "real attention."
- **T1.3 / F3 (MED) — Tighter softmax bounds.** `Softmax.compute_softmax_bounds` is sound but loose
  (interval box → many UNKNOWN). Implement `compute_softmax_bounds_tight` (Wei 2023 linear; currently a
  sampling placeholder). Tightening is sound-preserving — verify it still MC-contains **and** strictly
  narrows the box (fewer UNKNOWNs at fixed ε).

---

## Tier 2 — Scale and breadth

- **T2.1 — VNN-COMP 2025 full sweep on real compute.** Most benchmarks currently time out at the ~10 s
  smoke budget; they need **≥120 s per instance / cloud fan-out** to yield real verdicts. Use the now
  cache-accelerated GitHub Actions matrix (or Vanderbilt ACCRE SLURM array jobs / a cloud batch backend)
  to run the embarrassingly-parallel instance set. See `CLOUD_COMPUTE_OPTIONS.md` for the licensing +
  parallelization options (campus license / MLM token, `matlab-actions/setup-matlab`, Docker
  `mathworks/matlab`, ACCRE). Goal: minutes-not-hours iteration on the benchmark sweep.
- **T2.2 / F4 (LOW-MED) — VNN-COMP breadth (the 3 blocked categories).** `cctsdb_yolo`,
  `collins_aerospace`, `vit_2023` are fail-loud today and need new op semantics to import+verify:
  `Slice`/`Expand`/`Where`/`ScatterND`/`ArgMax`/large-`Split`/`Pow` (YOLO post-proc + Collins), and the
  **sound multi-token attention bound (T1.1)** for `vit_2023`. Add per-benchmark soundness tests as each
  lands.
- **T2.3 — Importer dynamic-layer soundness (continue).** Tests 27–30 cover ElementwiseProduct/Division/
  DynamicMatmul MC-containment; keep the **interval-division-across-zero** trap audited (bound-or-reject
  denominators straddling 0), and keep the rule that **unknown ops error at import** rather than loading a
  partial/identity net.

---

## Tier 3 — Adjacent workstreams & tech debt

- **T3.1 — SLM (sentiment / small-language-model) workstream.** A substantial separate effort (the
  `TODO_SLM-pass01..15` history) is **currently excluded** from PR #290 (its example/test files are
  gitignored, its CI bookkeeping references them as "former"). Decide its fate: either **complete it as its
  own follow-up PR** — with self-contained pass/fail soundness tests and a real CI wiring — or formally
  retire it. The orphan `tests/regression/run_regression_tests.m` (references the ignored
  `test_SLM_layer_outputs`) should be wired up or removed as part of that decision.
- **T3.2 / F5 (LOW) — Convert/retire the demonstration soundness tests.**
  `tests/soundness/{test_SLM_layers_soundness, test_SoftmaxLayer_reach_soundness}.m` *record* unsoundness
  rather than pass/fail and break under `runtests` (cross-`%%` scoping). Convert to self-contained
  pass/fail or drop them; their CI role is already covered by `test_transformer_soundness_regression.m`.
- **T3.3 — Pre-existing, out-of-scope subsystem bugs (flagged for maintainer).** These predate the PR and
  live in separate subsystems; each warrants its own issue:
  - `SignLayer.reach` is **unsound** (returns +1 only) — excluded from the exact-star whitelist here, but
    the method itself should be fixed or guarded.
  - `LinearNNCS`/`DLinearNNCS` `.verify` method-string gate + `.falsify` discards `U.contains`.
  - `verify_safety` parfor references an undefined `method`.
- **T3.4 — Formal soundness (stretch).** All current guarantees are **empirical (MC-containment)**. A
  rigorous follow-up could give formal soundness proofs for the value-hull multi-token bound (T1.1) and the
  softmax linear bounds (T1.3) — the priority targets for a paper-grade result.

---

## Suggested sequencing

1. **Land PR #290** (done once CI green on the head commit).
2. **Tier 0** in the first post-merge pass — especially **T0.2 (Node 24, deadline 2026-06-16)** and
   **T0.1 (branch protection)**; T0.3/T0.4/T0.5 are quick.
3. **T1.1 → T1.2** (sound multi-token → real ViT) as the flagship next feature, with T1.3 folding in to cut
   UNKNOWNs.
4. **T2.1 (cloud sweep)** in parallel with Tier 1 — it's infra, not engine, and unblocks honest VNN-COMP
   numbers; then **T2.2** breadth as op semantics land.
5. **Tier 3** opportunistically: the SLM decision (T3.1/T3.2) and the maintainer-owned pre-existing bugs
   (T3.3) as separate, well-scoped PRs/issues; T3.4 when a rigorous write-up is the goal.

**Verification habit (per MCP session):** `cd code/nnv; startup_nnv; NNVVERSION()` (expect `NNV v3.0.0`),
then `runtests('tests/soundness','IncludeSubfolders',true)` and
`runtests('tests/vnncomp25/test_vnncomp25_regression.m')` must stay green before and after each change.
