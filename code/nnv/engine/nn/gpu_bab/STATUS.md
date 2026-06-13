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

**Sound-or-unknown by construction:** `robust` only when every leaf's CROWN margins
exceed an FP slack; `unsafe` only from a concretely-evaluated misclassifying input;
`unknown` on budget. Never a wrong verdict (VNN-COMP −150 avoided).

## Known limitation (measured)

Input-split BaB does **not** scale to high-dim inputs: on MNIST (784-dim) it exhausts
the box budget → `unknown` for borderline eps. Input-split is efficient only for
low-dim inputs (acasxu, 5-dim). This is expected and matches the field.

## Roadmap (priority order)

1. **α-CROWN** — optimise the unstable-ReLU lower-slope `al` (currently fixed min-area)
   by projected gradient on the margin. Matches the LP triangle-relaxation optimum
   **without any LP**, and is *tighter* when intermediate bounds are jointly optimised.
   This is what makes CROWN competitive on MNIST; LP-free, GPU-friendly.
2. **β-CROWN ReLU-split BaB** — branch on unstable *neurons* (not input dims), fold the
   split constraint into propagation via an optimisable Lagrangian `β`. The scalable,
   complete method for MNIST-scale; sub-domains batch along the GPU dimension. Replaces
   input-split for high-dim inputs.
3. **gpuArray residency** — `gpu_bab_compile` to upload weights once as resident
   gpuArray; keep the BaB frontier on-device; gather only verdicts. (Needs a GPU
   instance to measure; license moves via ENI `eni-0a4349f2e1fa10e5f` to a g4dn.)
4. **Conv/ImageStar coverage** + dispatcher routing with Star fallback.

## Separate CPU win (from the LP-size measurement, not GPU-BaB)

cifar100 resnet solves **14,612 LPs of nVar=3072 (=input dim) but nCon=1** — i.e.
box-only LPs whose optimum is the closed-form corner `sum over i of max(c_i·lb_i, c_i·ub_i)`.
NNV calls `linprog` for these; replacing box-only predicates with the closed-form range
is a sound CPU speedup (LP-avoidance). Independent of the GPU engine.
