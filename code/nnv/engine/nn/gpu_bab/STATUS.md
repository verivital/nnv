# GPU-BaB engine — status & roadmap

LP-free, GPU-native bound-propagation + branch-and-bound verifier for NNV (the
research-backed alternative to GPU-accelerating NNV's LP-bound Star path, which is
impossible — see `../../../../../share_repo/research/gpu_polyhedral_reachability_research.md`).
Pure MATLAB, configurable precision (`'single'` default / `'double'`), gpuArray-ready.

## Built & CPU-validated (sound)

| file | what | validation |
|---|---|---|
| `nn_to_ops.m` | NN → affine/relu op list | op-eval == NNV evaluate (1e-15) |
| `gpu_bab_ibp.m` | sound interval bound prop | 0 violations / 60k samples |
| `gpu_bab_crown.m` | CROWN-IBP backward output bounds | 0 viol / 90k; 0.59× IBP; **sound on real mnistfc 256×4 (0/200k)** |
| `gpu_bab_crown_spec.m` | **batched** CROWN lower bound on a spec `C·f(x)` (pagemtimes over sub-domains) | sound (margin ≤ true min); **batched == single-box exactly (0.0)** |
| `gpu_bab_verify_robust.m` | single-pass robustness (robust\|unknown) | correct |
| `gpu_bab_bab.m` | input-split branch-and-bound | **sound `robust` (MC-confirmed) + sound `unsafe` (real counterexample)** |
| `gpu_bab_crown_alpha.m` | **α-CROWN** — optimise unstable lower slopes (dlgradient, normalised step + keep-best) | sound + **tightens (gain +8 to +35 on a 4-ReLU net)** |
| `gpu_bab_relu_split.m` | **β-CROWN-style ReLU-split BaB** — branch on neurons, clamp pre-activation bounds (sound; scales to high-dim) | sound `robust` (MC-confirmed); framework correct, see refinements |

**Sound-or-unknown by construction:** `robust` only when every leaf's CROWN margins
exceed an FP slack; `unsafe` only from a concretely-evaluated misclassifying input;
`unknown` on budget. Never a wrong verdict (VNN-COMP −150 avoided).

## BREAKTHROUGH (2026-06-14) — CROWN intermediate bounds: the engine now verifies

`gpu_bab_crown_tight` computes each layer's pre-activation bounds via a backward CROWN
pass (reusing earlier layers' tight bounds) instead of IBP. On the trained MNIST FC net
(2,304 ReLU) it tightens the robustness-margin lower bound by **~4 orders of magnitude**:

| eps | IBP-intermediate | CROWN-intermediate | engine verdict |
|---|---|---|---|
| 0.005 | −40,278 | **+5.84** | **robust** (root, 0.1s, MC-sound) |
| 0.01 | −84,714 | −270 | unknown (close) |
| 0.02 | −178,593 | −4,712 | unknown |

So the engine went from *useless* to a **sound verifier that proves real trained-MNIST
robustness** at small eps. Wired into the ReLU-split BaB (`intermediate='tight'`) +
BaBSR-lite gap branching (`i_pick_split_gap`).

### Remaining gap to competition parity (honest)

eps=0.01–0.05 on a 2,304-neuron net still returns `unknown` even with 800 BaB nodes — to
certify you must split away the slack of *hundreds* of unstable neurons.

**Two bounded experiments (2026-06-14) settled WHY it's stuck:**
- **Budget vs bound (3,000 nodes):** eps=0.01 reached **maxDepth=247** (fixed 247 neurons
  deep) and *still* `unknown`. So it is a **bound problem, not budget** — more nodes / GPU
  batching alone will NOT close it.
- **α-on-tight (final-pass, `useTight`):** buys only **~3%** (eps=0.01: −270 → −261). So the
  looseness is **not** in the final-layer slopes — it's in the **intermediate relaxations**.

**Therefore the only remaining lever is the FULL joint α-optimization** — α optimized *inside
every intermediate backward pass simultaneously* (the auto_LiRPA core), then GPU-batched BaB
for speed. That is essentially reimplementing α,β-CROWN in MATLAB: a **multi-week** project.

Current engine = a working, sound prototype that proves small-eps robustness. Competition
parity = the full joint α-CROWN + batched BaB. The box-LP fast-path (PR #355, merged) is the
immediate VNN-COMP win; the research independently concluded NNV's gap is coverage/method,
not compute.

## Known limitation (measured)

Input-split BaB does **not** scale to high-dim inputs: on MNIST (784-dim) it exhausts
the box budget → `unknown` for borderline eps. Input-split is efficient only for
low-dim inputs (acasxu, 5-dim). This is expected and matches the field.

### DECISIVE diagnostic (2026-06-14) — the bound is catastrophically loose on real nets

On a trained MNIST FC net (784→1024→512→256³→10, 2,304 ReLU), a confident digit
(clean margin **+15**), the CROWN robustness-margin lower bound is:

| eps | fixed-CROWN min | α-CROWN (60 it) min | proves at root |
|---|---|---|---|
| 0.005 | **−40,277** | −39,068 | no |
| 0.01 | −84,714 | −82,119 | no |
| 0.02 | −178,593 | −173,159 | no |

The bound is **~4 orders of magnitude** too loose (−40k vs the true +15), and 60
α-iterations move it only ~3%. **Root cause: this engine uses IBP for intermediate-layer
bounds**, which explode over a 6-layer net — and α-optimizing only the *final* spec
cannot fix loose *intermediate* bounds. No amount of α-iterations / BaB budget / branching
closes a 40,000× gap.

**Implication for the roadmap:** the missing piece is **CROWN intermediate bounds** — a
backward CROWN pass to bound *every* intermediate layer (the O(L²) core of real α-CROWN /
auto_LiRPA), with α optimized jointly across all those passes. That is a **major rework**,
not a tune. Until it lands, this engine cannot verify real-scale MNIST nets (it is sound —
returns `unknown` — never wrong). This confirms the research's conclusion that a
from-scratch competitive LiRPA verifier in MATLAB is a long investment, and that the
near-term VNN-COMP value is the **box-LP fast-path** (a real CPU win) + coverage/method
levers, not GPU compute.

## Roadmap (priority order)

DONE: **α-CROWN** (`gpu_bab_crown_alpha`) and the **β-CROWN ReLU-split framework**
(`gpu_bab_relu_split`) — both sound. The framework branches on neurons (scales past
input-split), but on a *hard random* MNIST net it still returns `unknown` within budget.
Refinements to make it *powerful* (priority order):

1. **α-CROWN bounds inside ReLU-split** — bound each BaB node with α-CROWN (clamped-bound
   variant) instead of fixed-slope CROWN. Tighter node bounds → far fewer splits → proves
   more without more budget. (The single biggest lever; α-CROWN already exists, needs the
   clamp plumbed through it.)
2. **Better branching (BaBSR-style)** — split the neuron whose split most improves the
   bound, not the first unstable one. Fewer nodes to a proof.
3. **Per-node counterexample search** + early termination — `gpu_bab_relu_split` currently
   searches counterexamples only once up front (so the huge-box test returned `unknown`
   instead of `unsafe`); search per node to falsify non-robust queries.
4. **Test on a TRAINED MNIST net** — the random net is a worst case (≈half the neurons
   unstable); a trained net has far fewer unstable neurons → many fewer splits. Likely the
   quickest way to demonstrate real `robust` verdicts.
5. **Batch the BaB frontier** through `gpu_bab_crown_spec` (per-node clamps as columns) —
   the GPU-parallel win once correctness is solid.
6. **gpuArray residency** — `gpu_bab_compile` to upload weights once as resident
   gpuArray; keep the BaB frontier on-device; gather only verdicts. (Needs a GPU
   instance to measure; license moves via ENI `eni-0a4349f2e1fa10e5f` to a g4dn.)
4. **Conv/ImageStar coverage** + dispatcher routing with Star fallback.

## Separate CPU win (from the LP-size measurement, not GPU-BaB)

cifar100 resnet solves **14,612 LPs of nVar=3072 (=input dim) but nCon=1** — i.e.
box-only LPs whose optimum is the closed-form corner `sum over i of max(c_i·lb_i, c_i·ub_i)`.
NNV calls `linprog` for these; replacing box-only predicates with the closed-form range
is a sound CPU speedup (LP-avoidance). Independent of the GPU engine.
