# Plan: Fix flaky soundness tests in `tests/soundness/`

Status: PLAN ONLY (no production code changed). Read-only analysis of:
- `code/nnv/tests/soundness/soundness_test_utils.m`
- `code/nnv/tests/soundness/test_soundness_MaxPooling2DLayer.m`
- `code/nnv/tests/soundness/test_soundness_ReluLayer.m`, `test_layer_soundness_harness.m`
- `code/nnv/engine/nn/layers/MaxPooling2DLayer.m`
- `code/nnv/engine/set/ImageStar.m`, `code/nnv/engine/set/Star.m`

Symptom: 1-4 of the 54 soundness tests red CI per run, non-deterministically; the canonical example is
`test_soundness_MaxPooling2DLayer/Test3_Multi_channelMaxPool` failing with
`Corner N not contained in output`. Same test passes most runs. `ci_allowed_failures.txt:24`
forbids allow-listing any soundness test, so these MUST be fixed.

---

## 1. Root-cause analysis: CHECK artifact, not a REACH-soundness gap

**Verdict: the flake is a containment-CHECK artifact, not a real reach-soundness bug.**
Two independent arguments, both backed by code lines:

### 1a. The math says a violation here is impossible if reach is correct
- `approx-star` for MaxPool is a *sound over-approximation* by construction. In
  `MaxPooling2DLayer.reach_star_approx` (`MaxPooling2DLayer.m:778-914`), when a pooling region has
  multiple max candidates it introduces a fresh predicate `y` with constraints
  `x_i - y <= 0` for every pool element `i` (`MaxPooling2DLayer.m:887-890`) and `y <= ub`
  (`MaxPooling2DLayer.m:878-882`). So `y >= max_i(x_i)` and `y <= ub`. The concrete output is
  `y* = max_i(x_i*)`, which is feasible (`y = max_i x_i` satisfies all `x_i - y <= 0` and, since
  `ub` is a sound per-region upper bound from `get_localBound`, `y <= ub`). **A correctly-computed
  approx-star output therefore always *contains* the concrete output.** A genuine "corner not
  contained" for approx-star or exact-star would be a real bug — but the evidence below shows the
  failures are at float-noise magnitude, i.e. the CHECK, not the SET, is wrong.

- For exact-star (Test1, Test4) the set is *exact*; the concrete output is a vertex of one of the
  output stars. A vertex lives on the boundary of `C*alpha <= d`, which is exactly where a
  tolerance-sensitive check misfires.

### 1b. The check itself is numerically fragile *and* under-specified
`verify_star_containment` / `verify_imagestar_containment` (`soundness_test_utils.m:8-179`) have a
"standard path" (`:61-82` and `:157-178`) taken whenever `n_pred <= n_out` (the usual case for these
small outputs):

```
alpha_solve = V_basis \ residual;      % or pinv if rank-deficient   (line 64-66 / 159-163)
reconstruction_error = norm(V_basis*alpha_solve - residual);          (line 70 / 166)
if reconstruction_error > tol -> contained = false                    (line 71 / 167)
contained = all(C*alpha_solve <= d + tol)                             (line 80 / 176)
```

Three concrete defects, each able to flip a true-contained point to `false`:

1. **The box bounds `predicate_lb <= alpha <= predicate_ub` are never enforced in the standard
   path.** Only the LP path (`:42-60`, `:124-156`, taken when `n_pred > n_out`) passes `lb/ub` to
   `linprog`. The MaxPool relaxation predicate `y` is bounded by `pred_lb=lb_region`,
   `pred_ub=ub_region` (`MaxPooling2DLayer.m:880-881`), and its *only* `C` rows are the relaxation
   inequalities. The least-squares `V_basis \ residual` minimizes reconstruction norm with **no
   knowledge of those bounds**, so for a boundary corner it can return an `alpha` that reconstructs
   the point but sits slightly outside the true predicate box — and then `C*alpha <= d + tol` is evaluated at
   the wrong `alpha`. This is *exactly* the geometry of a corner sample.

2. **Boundary corners make `C*alpha <= d` marginal at the `1e-6` level.** The active relaxation
   constraint is `x_i* - y = 0` at the corner. `V_basis \ residual` solves an over-determined LS;
   its residual is `O(cond(V_basis) * eps * ||residual||)`, typically `1e-9`..`1e-5` depending on
   the random `V` (Test3 builds the center with `rand(4,4)`, `test_soundness_MaxPooling2DLayer.m:51-52`).
   That residual feeds straight into `C*alpha`, so a constraint that should read `<= d` reads
   `<= d + (1e-9..1e-5)`. With a *fixed* `tol = 1e-5` inside `verify_layer_soundness_*`
   (`soundness_test_utils.m:316,412`) but a **default `tol = 1e-6` inside the containment helpers**
   (`:20,98`) — and the helpers are called *with* the passed `tol`, so 1e-5 here — the margin is
   thin and occasionally exceeded. The single tight absolute `tol` does not scale with the
   magnitude of `residual` or `d`.

3. **Non-determinism comes from unseeded `rand`.** Test3 center is `rand(4,4)` per channel
   (`test_soundness_MaxPooling2DLayer.m:51-52`); the random samples in `sample_imagestar`
   (`soundness_test_utils.m:239`) are unseeded. The random center changes *which* pooling regions
   tie (and thus the output set's constraint geometry and `cond(V_basis)`), and the unseeded sampler
   changes which points land near a face. So the same test passes on most seeds and fails on the few
   where a corner/sample lands on a face with an ill-conditioned `V_basis`. The corner enumeration
   itself (`get_predicate_corners`, `:259-295`) is deterministic — the flake is the *data*
   (center + samples), not the corner indexing.

### Why the naive `lsqlin` swap made it WORSE (the reverted attempt)
Replacing `V_basis \ residual` with constrained `lsqlin` (box + inequality constraints, equality
`V*alpha = residual`) failed *deterministically even for exact-star points* because:
- `lsqlin` defaults to the `'interior-point-convex'` (or `'trust-region-reflective'`) algorithm with
  *relatively loose* default `OptimalityTolerance`/`ConstraintTolerance` (`~1e-8..1e-6`) and stops at
  an interior-ish iterate; for a point that is a *vertex* (exact-star) the true `alpha` is at a
  corner of the feasible box, where interior-point methods converge slowly and return an `alpha`
  whose reconstruction error exceeds the tight `1e-6`..`1e-5` check — so `reconstruction_error > tol`
  fires and the point is rejected. In other words, swapping the *solver* without *loosening and
  band-structuring the acceptance test* trades "occasionally too strict on boundary noise" for
  "reliably too strict on every vertex." The fix is not a different solver alone; it is a feasibility
  test with a **tolerance band that matches the competition's notion of equality** (Section 3) plus
  the **box bounds** (Section 4).

### What measurement confirms this (for the main agent to run)
Instrument `verify_imagestar_containment` to, on a `false` return, print
`reconstruction_error`, `max(C*alpha - d)`, and `max(alpha - pred_ub, pred_lb - alpha)`.
**Prediction:** on the flaky failures these are all `< 1e-4` (float/LP noise), never `O(0.1)`
(which is what a *real* unsound reach looks like — cf. the LayerNorm `~0.875` overflow noted in
`test_layer_soundness_harness.m:9-11`). If any violation is `>= 1e-3`, escalate to a real reach bug
and stop applying the tolerance fix.

---

## 2. Determinism vs fuzzing — recommendation

**Recommendation: deterministic fuzzing — seed `rng` to a small fixed set of seeds and loop.**

Keep the random *coverage* (it is what catches real unsoundness, per the harness rationale at
`test_layer_soundness_harness.m:7-14`) but make every CI run reproducible.

| Option | Pro | Con |
|---|---|---|
| (A) Leave unseeded `rand` | Maximal long-run coverage; different inputs each CI run | Flaky/non-reproducible; a failure can't be reproduced locally; reds PRs unrelated to the change |
| (B) Single `rng(seed)` per test top | Fully reproducible; trivial 1-line change | Only ONE input per test forever — coverage collapses to a single point; a latent bug on other seeds is invisible |
| (C) **Loop over K fixed seeds (deterministic fuzzing)** | Reproducible AND broad: K inputs per test, all replayable; a failure prints the exact seed | Slightly more runtime (~Kx the sampling, not the reach); needs a tiny loop refactor |

Choose **(C)**. Concretely: at the top of each test (or in a shared helper), iterate
`for seed = SOUNDNESS_SEEDS` with `SOUNDNESS_SEEDS = [0 1 2 7 42]` (5 seeds), call `rng(seed)`, build
inputs, and assert containment; on failure include the seed in the message
(`'... failed (seed=%d): %s'`). This gives reproducibility (CI is now a pure function of the code)
while preserving fuzz breadth. The existing harness already proves the pattern works: it uses a
single `rng(7)` (`test_layer_soundness_harness.m:16`) and is *not* flaky — generalize that to a seed
loop so we don't lose coverage to a single point.

Guardrail for the constraint in `ci_allowed_failures.txt`: deterministic seeding does **not** weaken
the gate — an unsound reach fails on *every* seed, so a fixed seed set still catches structural
unsoundness; it only removes the float-noise flake.

---

## 3. Tolerance — grounded in the VNN-COMP 2026 rules

VNN-COMP accepts a counterexample witness when onnxruntime replay matches NNV's outputs to
**< 1e-3 relative error** and constraint satisfaction to **1e-4 absolute error**. The competition
itself does not treat two numbers closer than that as "different." A containment test that rejects a
point because a reconstruction/LP residual is `2e-6` is therefore **stricter than the standard the
verifier is judged by** — it is testing the LP solver's float noise, not soundness.

**Recommendation: replace the single absolute `tol = 1e-6/1e-5` with a paired absolute+relative
band, matched to the competition:**

```
abstol = 1e-4;     % matches VNN-COMP constraint absolute tolerance
reltol = 1e-3;     % matches VNN-COMP output relative tolerance
```

Acceptance becomes, for each scalar comparison value `v`:
`ok(v_lhs, v_rhs) := (v_lhs <= v_rhs + abstol + reltol*abs(v_rhs))`,
applied to (i) reconstruction `|V*alpha - residual| <= abstol + reltol*|residual|` componentwise and
(ii) constraints `C*alpha <= d + abstol + reltol*abs(d)` and box bounds with the same band.

Why this catches real unsoundness while ignoring noise:
- A *real* unsound reach (a missing reachable output / an ONNX-vs-NNV inference difference) is
  **structural and large** — the LayerNorm/SiLU bugs that motivated the harness overflowed by
  `~0.1`..`0.875` (`test_layer_soundness_harness.m:9-11`), `>> 1e-3`. The band leaves a 3-order-of-
  magnitude gap between "float noise" and "real bug," so detection power is preserved.
- A boundary corner whose LS/LP residual is `1e-6`..`1e-5` is now *inside* the band and passes,
  killing the flake.

Risk both ways (state explicitly):
- **Too tight** (status quo `1e-6`): spurious failures from LP/float noise — the current bug.
- **Too loose** (e.g. `1e-2`): could mask a *small but real* unsoundness or a genuine inference
  difference below `1e-2`. Mitigation: stay at the competition's own numbers (`1e-4`/`1e-3`); do not
  go looser without evidence. Keep the harness's interval check at its existing `1e-6`
  (`test_layer_soundness_harness.m:17`) for the layers it covers if desired, OR unify on `1e-4/1e-3`
  — but if unifying, first run the measurement in Section 1 to confirm no covered layer relies on the
  tighter bound to catch a real bug.

---

## 4. The containment-check correctness fix

### Prefer the existing, exact NNV method where available
NNV already ships exact containment via `Polyhedron` (MPT) feasibility:
- `Star.contains(s)` (`Star.m:1119-1143`): builds `Polyhedron('A',C,'b',d,'Ae',V(:,2:end),
  'be', s - V(:,1),'lb',predicate_lb,'ub',predicate_ub)` and returns `~isEmptySet`. **This already
  enforces the predicate box bounds and the equality** — i.e. it does exactly what the hand-rolled
  standard path forgets.
- `ImageStar.contains(image)` (`ImageStar.m:611-637`): reshapes and delegates to `Star.contains` via
  `toStar`.

These are the *correct* feasibility formulation. The hand-rolled helper exists because (a)
`Polyhedron`/`isEmptySet` has **no tolerance band** (a boundary vertex can be reported empty under
MPT's own internal tolerance, reintroducing the same flake — note `Star.isEmptySet` at
`Star.m:1146-1162` was itself written to work around an MPT `isEmptySet` bug), and (b) it adds an MPT
dependency in the hot path. So we do not blindly swap to `Star.contains`; we fix the helper to be the
*same feasibility test* but with the competition tolerance band and the box bounds, using `linprog`
(already a hard NNV dependency via `lpsolver`) rather than MPT.

### The robust formulation (one path for all cases)
Replace the two-branch (`n_pred > n_out` LP vs else LS) logic with **a single LP feasibility test
that always includes the box bounds and a tolerance band**, written as an `lpsolver`/`linprog`
feasibility problem. For a point `p`, center `c`, basis `V`, constraints `C*alpha <= d`, box
`[plb, pub]`:

Find `alpha` s.t.
```
plb <= alpha <= pub                                  (box — the part the standard path dropped)
C*alpha <= d + (abstol + reltol*abs(d))              (banded inequalities)
-(abstol + reltol*abs(r)) <= V*alpha - r <= (abstol + reltol*abs(r))   r = p - c  (banded equality)
```
Implement the banded equality as two inequality blocks (so the whole thing is a single `linprog`
feasibility LP with `f = 0`, no `Aeq`): this avoids `linprog` requiring an *exact* `Aeq*alpha = beq`,
which is what made both the original equality-LP path and the `lsqlin` attempt brittle on vertices.
`contained = (exitflag == 1)`. Using the slack-equality form (inequalities only) means a vertex point
is *interior* to the banded feasible region, so the LP converges cleanly regardless of solver
algorithm — directly fixing the `lsqlin`/`Aeq` brittleness from Section 1.

This is minimal and low-risk: same variables, same `linprog` already used at
`soundness_test_utils.m:52-53`, just (1) always include `lb/ub`, (2) drop the exact `Aeq` in favor
of two banded inequality blocks, (3) use the abstol+reltol band everywhere a `< tol` comparison
exists today.

Keep a fast pre-check: if `n_pred == 0`, the existing center-equality test
(`soundness_test_utils.m:30-34, 103-107`) stays, but with the band:
`norm(point - c) <= abstol + reltol*norm(c)`.

---

## 5. Concrete step-by-step implementation plan (least-risky first)

> All steps confined to `code/nnv/tests/soundness/` (test code). **No `engine/` change** — the
> reach is sound; only the test harness is wrong.

**Step 0 — Measure (confirms diagnosis, gates the rest).**
Temporarily add diagnostic prints to `verify_imagestar_containment`'s `false` branches
(`reconstruction_error`, `max(C*alpha-d)`, box overflow). Run
`test_soundness_MaxPooling2DLayer` in a loop (e.g. 200 iterations, unseeded) until a failure occurs;
record the magnitudes. **Expected: all `< 1e-4`.** If any `>= 1e-3`, STOP and open a reach-soundness
investigation instead. (Revert the prints after.) — *risk: none (read-only instrumentation).*

**Step 1 — Add abstol/reltol band to the containment helpers (no signature change).**
In `soundness_test_utils.m`, introduce `abstol = 1e-4; reltol = 1e-3;` (derive from the passed `tol`
for backward-compat: `abstol = max(tol, 1e-4)`), and replace every `< tol` / `<= d + tol` /
`norm(...) < tol` comparison in `verify_star_containment` (`:30-82`) and
`verify_imagestar_containment` (`:103-178`) with the banded form from Section 3/4. — *risk: low; only
loosens acceptance, monotonic, cannot introduce a false PASS larger than the band.*

**Step 2 — Make the standard path enforce box bounds via the single-LP feasibility (Section 4).**
Collapse the `if n_pred > n_out ... else ...` branches into the one banded `linprog` feasibility LP
(inequality-only, box included). This removes the bounds-ignoring `V\residual` path entirely. — *risk:
low-moderate; new code path, but it strictly subsumes the old LP path which already worked.*

**Step 3 — Deterministic fuzzing seeds.**
Add `SOUNDNESS_SEEDS = [0 1 2 7 42]` and wrap each test's input construction + assertion in a
`for seed = SOUNDNESS_SEEDS; rng(seed); ...; end` loop, propagating `seed` into failure messages.
Apply to all 44 tests that use unseeded `rand` (start with `test_soundness_MaxPooling2DLayer.m`, then
the rest). A shared helper `soundness_test_utils.seeds()` keeps it DRY. — *risk: low; pure test
harness; increases coverage.*

**Step 4 — Optional: prefer `Star.contains`/`ImageStar.contains` as a cross-check (not the primary).**
For exact-star tests, optionally also assert via `ImageStar.contains` to triangulate the new helper
against NNV's own (MPT) method on the cases where MPT is available. Skip if MPT not installed. — *risk:
low; additive.*

### Validation plan (confirm 0 flakes without masking a real bug)
1. **Stress loop:** run `runtests('test_soundness_MaxPooling2DLayer')` (and the full
   `tests/soundness` suite) **200x unseeded** *before* Step 3 with the new banded helper — expect
   0 failures. This proves the *check* fix alone removes the flake even under randomness.
2. **Determinism:** after Step 3, run the suite 5x — byte-identical pass/fail (it is now a pure
   function of code).
3. **Mutation / no-masking test (critical):** inject a deliberately UNSOUND reach (e.g. shrink the
   MaxPool relaxation `ub` by `0.1`, or in a scratch copy widen a sample beyond the set by `0.1`) and
   confirm the banded check **still fails** — proving the `1e-4/1e-3` band did not blind the gate to a
   real `O(0.1)` violation. Revert.
4. **Harness parity:** run `test_layer_soundness_harness.m` (already seeded, already banded at
   interval level) to ensure unification of tolerances (if done) doesn't regress its detection of the
   LayerNorm/SiLU-class bugs.

### Risks & mitigations
- *Risk:* the band hides a small (`~1e-3`) real inference difference. *Mitigation:* anchor exactly to
  the competition numbers (no looser); Step 0 measurement + Step 3 mutation test bound the blind spot.
- *Risk:* `linprog`-only feasibility differs subtly from MPT `Polyhedron`. *Mitigation:* Step 4
  cross-check on exact-star where MPT is present.
- *Risk:* seed loop increases CI time. *Mitigation:* 5 seeds, and the cost is in *sampling*, not
  *reach*; cap `n_samples` if needed (the gate is corners + a modest sample count, not exhaustive).
- *Risk:* touching 44 files. *Mitigation:* land Steps 1-2 (the actual fix) first to stop the bleeding;
  Step 3 (seeds) can follow per-file. The check fix alone makes the suite robust even while still
  unseeded (validation item 1).

---

## TL;DR
The flake is a **containment-CHECK artifact**, not unsound reach: the helper's standard path uses
`V\residual` least-squares that (a) **ignores the predicate box bounds** and (b) evaluates
`C*alpha <= d + 1e-6` on boundary corners where LP/float residual (`1e-9..1e-5`) tips a marginal
constraint over a too-tight tolerance; unseeded `rand` makes which-corner-tips non-deterministic. Fix:
(1) **deterministic fuzzing** over a fixed seed set, (2) **abstol 1e-4 / reltol 1e-3** band matching
VNN-COMP's own equality tolerances, (3) a **single banded `linprog` feasibility LP that includes the
box bounds** (mirroring `Star.contains` but with a tolerance band and no exact `Aeq`, which is why the
prior `lsqlin` swap failed). All changes are test-only; no `engine/` reach change. Validate with a
200x stress loop and a mutation test that confirms an injected `O(0.1)` unsoundness still fails.
