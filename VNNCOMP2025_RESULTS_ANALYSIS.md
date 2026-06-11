# VNN-COMP 2025 results — competitive analysis for NNV (where we lose, where to gain)

*Hard numbers from the official results repo (github.com/VNN-COMP/vnncomp2025_results,
SCORING-SMALL-TOL/latex), 2026-06-11. This is the empirical basis for
[`VNNCOMP2026_STRATEGY.md`](VNNCOMP2026_STRATEGY.md). The narrative report is arXiv:2512.19007.*

## Final 2025 leaderboard (overall score) — NNV was 6th of 7

| # | Tool | Score | Method (short) |
|---|------|------:|----|
| 1 | α-β-CROWN | **1600.0** | GPU linear-bound propagation (CROWN) + branch-and-bound (β-CROWN / GCP-CROWN) |
| 2 | NeuralSAT | 1360.1 | GPU; DPLL(T) / abstraction + clause learning |
| 3 | PyRAT | 1218.5 | abstract interpretation (zonotope/polytope), CPU |
| 4 | **CORA** | 954.6 | **MATLAB** zonotope/polynomial-zonotope reachability |
| 5 | nnenum | 740.8 | star/zonotope, CPU, ReLU case-splitting |
| 6 | **NNV** | **697.3** | **MATLAB** star/ImageStar/zono/abs-dom reachability |
| 7 | SobolBox | 529.0 | sampling/Sobol |

**The gap is not "we're MATLAB/CPU."** CORA (954.6) and nnenum (740.8) use the *same* family of
reachability tech as NNV and beat it. The gap is **SAT-finding, penalty leakage, and speed** — all fixable.

## Diagnostic tables (the why)

| Metric | α-β-CROWN | NeuralSAT | PyRAT | CORA | nnenum | **NNV** | SobolBox |
|---|---:|---:|---:|---:|---:|---:|---:|
| Benchmarks participated | 16 | 16 | 15 | 14 | 11 | **15** | 9 |
| Instances verified (total) | 2575 | 2425 | 2090 | 1854 | 1596 | **1082** | 1057 |
| **SAT found** | 1082 | 1055 | 1018 | 946 | 786 | **354** | 317 |
| UNSAT proved | 1493 | 1370 | 1072 | 908 | 810 | **728** | 740 |
| **Incorrect / missing-CE** | 0 | 17 | 2 | 20 | 0 | **19** | 437 |
| Startup overhead (s) | 6.0 | 5.9 | 6.0 | 10.9 | 0.9 | **14.2** | 2.7 |

## The four NNV-specific takeaways (ranked by score impact)

1. **🔴 SAT-finding is the #1 weakness — 354 vs ~1000.** NNV finds **3× fewer counterexamples** than the
   leaders. Its runner does only ~100–500 *random* samples (`falsify_single`/`create_random_examples`),
   while the field uses **gradient/PGD adversarial attacks**. SAT cases are CHEAP points (no reachability
   needed) — closing even half this gap is a large, low-cost score gain. **→ add a PGD/gradient
   falsifier** (autodiff through the `dlnetwork` toward the unsafe `HalfSpace`), multi-start, run it FIRST.
2. **🔴 Penalty leakage — 19 incorrect / missing-CE.** NNV "incorrect" is almost certainly **invalid or
   missing counterexample witnesses** (claimed SAT but the written witness didn't replay/validate), not
   wrong reachability — matching the runner's own note ("penalties last year… wrong counterexample
   writing?"). Each costs heavily. **→ validate every SAT witness by replaying it through the network
   (onnxruntime/`evaluate`) before emitting `sat`; fix the needReshape/flatten-order witness encoding** —
   never emit `sat` without a confirmed counterexample.
3. **🔴 Slowest startup (14.2 s).** `startup_nnv` + path/`javaaddpath`/dependency-check overhead is the
   worst in the field; on short-timeout instances this is pure lost budget. **→ a lean
   `prepare_instance.sh`/`run_instance.sh` path that does the heavy setup ONCE (install phase), caches a
   warmed MATLAB/NNV state, and trims per-instance startup.**
4. **🟠 UNSAT (proving) is closest to the pack (728) but still trails 2×.** Reachability proving is NNV's
   relative strength; the leaders win via **GPU bounds + branch-and-bound refinement** (prove the easy part
   with cheap bounds, split only where needed) rather than one expensive pass. **→ cheap-to-precise
   abstraction ladder + input/ReLU splitting refinement; explore GPU bound propagation.**

**Coverage (15/16 benchmarks) is already competitive — do NOT over-invest in breadth; invest in SAT-finding,
witness validity, and speed.**
