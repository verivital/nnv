# NNV — VNN-COMP 2026 Scoring-Maximization Strategy

*How NNV moves up from 6th/7 (697.3 pts, 2025). Grounded in the hard 2025 scoreboard
([`VNNCOMP2025_RESULTS_ANALYSIS.md`](VNNCOMP2025_RESULTS_ANALYSIS.md)), the 2026 rules, the 2025 report
(arXiv:2512.19007), and deep research on the winning tools' methods. Companion:
[`VNNCOMP2026_READINESS.md`](VNNCOMP2026_READINESS.md) (benchmarks/vnnlib2/deadline). **Deadline: tool
submission June 30, 2026 AoE.***

## The score math → what actually moves the needle

Per benchmark: **correct sat/unsat = +10, INCORRECT = −150, timeout/error/unknown = 0**, no time bonus.
Benchmark score = % of the max tool's score on that benchmark; overall = sum of per-benchmark %s; ties → total
runtime. So the levers, in order of leverage **for NNV given the 2025 data**:

1. **Find more SAT** — NNV got **354 SAT vs ~1000** for the leaders. SAT instances are *cheap points* (no
   reachability needed), and NNV's falsifier is the weakest in the field (random sampling only). **This is
   the single highest-ROI fix.**
2. **Stop leaking −150** — NNV had **19 incorrect/missing-CE** (invalid SAT witnesses). At a 16× penalty, a
   handful of these erase dozens of correct answers. **Validate every witness before emitting `sat`.**
3. **Be faster** — NNV had the **worst startup (14.2 s)** and times out on big nets; more instances solved
   within timeout = higher %, and runtime is the tiebreak.
4. **Prove more UNSAT** — NNV's 728 is its relative strength but still ~2× behind; close it with
   cheap-bounds + refinement (not one expensive exact-star pass).

**NNV's soundness gate is a competitive ASSET** — the −150 rule punishes the aggressive; "never certify
not-robust/unsafe from an over-approximation, return `unknown` instead" is exactly the right posture.

## The portfolio pipeline (per instance)

Replace the current "falsify-100-random → reach-method-list" with a **time-budgeted, falsify-first,
cheap-to-precise** pipeline, run with **falsification and reachability racing in parallel** (`parfeval`),
first sound verdict wins, the loser is killed:

```
budget T (per-instance timeout)
├─ FALSIFY worker  (target: SAT fast)         ┐ race; first {sat|unsat} wins,
│   FGSM warm-start → PGD (multi-start) →      │ kill the other. Always VALIDATE
│   CW fallback → random.  VALIDATE witness.    │ a SAT witness before emitting.
├─ PROVE worker    (target: UNSAT fast)        ┘
│   interval/zono → abs-dom → approx-star,  with
│   ABSTRACTION-REFINEMENT (split top-k unstable
│   neurons) instead of jumping to exact-star.
└─ at T−ε: emit `unknown` (0 pts) — NEVER a partial/unsound verdict.
```

Never run `exact-star` on large nets (it's the exponential trap that caused 2025 timeouts; the 2025 fix
already leads the compute-bound categories with `approx-zono`/`abs-dom`).

## Pillar 1 — Falsification (the #1 score lever: close 354 → ~900 SAT)

NNV's `falsify_single` does only ~100–500 *random* samples, no gradients. Add a real adversarial attack:
- **Loss = distance to the unsafe region:** for unsafe `U = {y : Hg·y ≤ g}`, minimize
  `max(0, margin-to-violation)` over `x ∈ [lb,ub]` via MATLAB autodiff (`dlarray`/`dlfeval`/`dlgradient`)
  through the `dlnetwork`.
- **FGSM warm-start → PGD refine:** 1-step FGSM on ~20 seeds → keep the top violating candidates → 50-step
  PGD (Adam, projected to `[lb,ub]`). **Multi-start in `parfor`** (each worker its own net copy).
- **Adaptive ε** by domain (images ≈ 0.1–0.5; flat features ≈ 0.01–0.1 of the range).
- **CW attack** only as a last resort on hard, highly-constrained instances (nn4sys, sat_relu).
- **Per-benchmark-class tuning** (`NNV_NRAND_FALSIFY` + PGD steps): small FC (ACAS/safenlp) — few samples,
  fast PGD; large images (cifar100/tinyimagenet/vgg) — many samples, more PGD steps.
- **Run falsification FIRST and in parallel with reach**, not just as a 0.5 s pre-pass.
- *Expected:* the literature + the 354-vs-1000 gap suggest PGD roughly 2–3×'s SAT yield — the biggest single
  point gain available, at modest effort.

## Pillar 2 — Validity & soundness (drive the 19 penalties → 0)

- **`validate_witness.m`:** before writing `sat`, **re-evaluate the candidate `x` through the network and
  confirm the property is violated with a 1e-6 input-constraint margin** (issue #2 requires zero-tolerance
  on input constraints; output replay tolerated to 1e-3 rel / 1e-4 abs). If it fails → emit **`unknown`,
  never `sat`.**
- **Fix witness encoding:** the `needReshape`/flatten ordering must match the ONNX input order (the runner
  already has `needReshape==3` handling — extend + test it for every benchmark; the 2025 penalties were
  "wrong CE writing").
- Keep the **`[42]` exact-star gate** (no not-robust from an over-approximation) — it's why NNV's incorrect
  count comes only from witnesses, not reachability.

## Pillar 3 — Imprecise-but-sufficient + refinement + GPU (close UNSAT, get fast)

- **CEGAR-style abstraction-refinement (`NN_reach_refinement.m`):** start with a cheap over-approx
  (interval/zono/abs-dom); if undecided at a time threshold, **split the top-k unstable ReLU neurons** (rank
  by activation range straddling 0), spawn the `2^k` phase-fixed sub-Stars in `parfor`, combine by
  disjunction (SAT if any child SAT, UNSAT if all children UNSAT). This is the middle ground NNV lacks
  between loose-zono and exponential exact-star.
- **Adaptive `relaxFactor`** per benchmark: 0.8–0.9 for large nets (ViT/VGG/YOLO), 0.0–0.2 for small (ACAS).
- **GPU bound propagation (high-impact, optional):** MATLAB `dlnetwork`/`dlarray` run on GPU. Use it for
  **PGD (Pillar 1) and CROWN-style linear bounds** on large CNNs — this is exactly where α-β-CROWN wins.
  Gate behind `NNV_USE_GPU`; CPU fallback. (Note: the p3 box has only 61 GB RAM — see Platform.)
- **GCP-CROWN-style cutting planes:** NNV's Star already supports extra constraints (`C·a ≤ d`); append
  tightening cuts from ReLU stability after a loose pass.
- **Bounds caching:** ACAS-Xu runs the *same net* over many input sets — cache per-layer interval bounds
  across a benchmark's instances.
- **Bandit method selection:** record `(benchmark, method, status, time)` during the pre-June-30 cloud
  sweep; reorder `reachOptionsList` per benchmark to put the historically-winning method first.

## Speed & parallelism (NNV was the slowest to start)

- **Keep `prepare_instance.sh` minimal** (it is); do heavy setup ONCE in `install_tool.sh`. Trim
  `startup_nnv` per-instance cost; consider a pre-warmed MATLAB/NNV state.
- **Use the 64 cores:** `parfor` over (a) multiple vnnlib specs (already partial), (b) falsification
  multi-start, (c) refinement sub-problems. Avoid nested parpools (the known hang).
- **Single-thread MATLAB where parpool isn't helping** (drop JIT/pool overhead on tiny instances).

## Platform & submission harness

- **Platform: `m5.16xlarge` CPU (64 vCPU, 256 GB)** as primary — NNV is CPU/LP/star-based and the large
  ViT/VGG/TinyImageNet nets need the RAM the GPU box (61 GB) lacks. **Re-evaluate GPU (`p3.2xlarge`) only if
  the PGD-on-GPU falsifier + CROWN bounds prove decisive** and fit in 61 GB.
- **Scripts:** `install_tool.sh` (deps/toolboxes/licenses + warm-up), `prepare_instance.sh` (≤10 min, no
  analysis), `run_instance.sh` (single-word result + validated SAT witness). **Start unofficial dry-runs on
  the TUM repeatability site NOW** to surface install/format errors early.

## VNN-LIB 2.0 (extended track) — use the Python packages (your call)

Per your direction, **bridge via the existing Python `vnnlib`/`vnnlib2` packages** rather than reimplement
the 2.0 grammar in MATLAB by June 30: extend the vendored `onnx2nnv_python` pattern with a
`vnnlib2nnv.py` that parses (pip `vnnlib` v1.0.2 handles 2.0 + `.vnnlib.gz`) and emits a 1.0-style
bounds/spec `.mat` NNV already consumes. Defer the deep `declare-hidden` (intermediate-tensor) and
multi-network (`isomorphic/monotonic_acasxu`) cases unless time allows — the **regular track (24, all VNN-LIB
1.0) needs no 2.0 parser at all** and is where most of the score is.

## Sequenced engineering plan (effort vs the June 30 deadline)

**🔴 Tier 0 — quick, highest-ROI (do first):**
1. **PGD/FGSM falsifier + multi-start `parfor`** (Pillar 1) — the 354→~900 SAT gap. *Biggest point gain.*
2. **`validate_witness.m`** + witness-encoding fixes (Pillar 2) — kill the 19 −150 leaks.
3. **Submission scripts + TUM dry-run** — surface install/format issues now.
4. **Full cloud sweep** (`vnncomp-sweep.yml`, official repo) at realistic timeout → the real regular-track
   scorecard + the bandit/method data.

**🟠 Tier 1 — medium:**
5. Per-instance **falsify∥reach race** (`parfeval`) + timeout-to-`unknown` guard.
6. Auto-tuned `relaxFactor` + bandit method ordering + bounds caching.
7. Startup-overhead reduction; broader `parfor` parallelism.

**🟡 Tier 2 — deeper (time-permitting):**
8. **Abstraction-refinement loop** (`NN_reach_refinement.m`, multi-neuron split) — the UNSAT closer.
9. **GPU** bound propagation + PGD on GPU.
10. **VNN-LIB 2.0 Python bridge** for the 4 extended-track 2.0 benchmarks.

**Bottom line:** the data says the win condition is **falsification + witness validity + speed**, not new
verification theory. Tier 0 alone (better SAT-finding and zero penalty leakage) should move NNV materially
up the table; the soundness gate keeps the −150 column at zero, which is where the aggressive tools bleed.
