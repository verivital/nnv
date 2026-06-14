# GPU-BaB engine — design constraints (layer coverage + GPU memory model)

Two architectural realities drive this engine's design. Both were raised as risks
and are addressed here so the implementation respects them from the start.

## 1. Why a separate engine at all (and the layer-coverage cost it implies)

NNV's existing Star/ImageStar `reach` path decides ReLU splits with a **per-neuron
LP** (`linprog`/GLPK). `engine/utils/lpsolver.m` *gathers any gpuArray back to the
host and forces double* before solving — i.e. **there is no GPU LP in MATLAB**. So
the existing path cannot be GPU-accelerated in place; GPU execution **requires** an
LP-free bound-propagation engine (IBP/CROWN: pure linear algebra, batched matmul).
This is the deliberate Option-B tradeoff from the plan.

The price of LP-free + GPU is that this engine does **not** reuse NNV's per-layer
`reach` methods — it needs its own bound-prop rule per layer type. That is the cost
flagged in review. Three things bound it:

- **The affine family is ONE rule, not N.** Conv2D, BatchNorm (eval mode), AvgPool,
  and FullyConnected are all *linear* operators. In CROWN/IBP they share a single
  backward rule — `A := A·op`, `d += A·bias` — differing only in the linear operator
  and its adjoint (transpose) for the backward pass. Conv's adjoint is a transposed
  convolution; AvgPool's is an upsample-and-scale; BatchNorm-eval is a per-channel
  diagonal scale+shift. So "many layer types" collapses to: the affine family (one
  rule, a few operators) + a small set of nonlinearities.
- **Only a handful of nonlinearities.** ReLU (done). MaxPool needs a relaxation
  (treatable like ReLU over its argmax window) — the main remaining nonlinear rule
  for image benchmarks. Flatten/Reshape is identity on values.
- **The Star path is the safety net → nothing is lost.** `nn_to_ops` ERRORS on any
  unsupported layer; the dispatcher catches that and falls back to the existing Star
  pathway. GPU-BaB is therefore **purely additive and opt-in per benchmark** — we
  grow coverage starting with the benchmarks that need it most (big conv image nets
  where Star is too slow), and every other benchmark keeps its current path. No
  benchmark regresses because of a coverage gap.

**GPU-target benchmark layer inventory (what we actually must cover, not all of NNV):**
Conv2D, BatchNorm, AvgPool/MaxPool, Flatten, FullyConnected, ReLU. Everything but
MaxPool is affine or reshape; ReLU is done. That is the bounded scope.

## 2. GPU memory + type model (the performance crux)

Confirmed against the R2026a Parallel Computing Toolbox docs ("Run MATLAB Functions
on a GPU", "Measure and Improve GPU Performance"):

- **Type must be explicit, and `single` is the default.** "Most GPUs perform
  calculations faster in single precision." Consumer-class GPUs (incl. **T4 / A10G**,
  our likely g4dn/g5 targets) have **24–64× more FP32 units than FP64**; datacenter
  A100/H100 only 2×. Query `gpuDevice().SingleDoubleRatio`. → **R1 maps directly:**
  `precision='single'` default (fast on T4/A10G); `'double'` is the opt-in tighten
  knob, expensive on consumer GPUs.
- **Host↔device transfer is THE bottleneck.** "Transferring data back to the CPU can
  be costly... limit the number of times you transfer data between host memory and
  the GPU." Design rules that follow:
  1. **Compile weights to the GPU ONCE.** A `gpu_bab_compile(net, precision)` step
     uploads every weight/bias as a resident `gpuArray` and returns an op list that
     already holds device arrays. The per-call `cast(op.W,...)` in the current
     single-box `gpu_bab_ibp`/`gpu_bab_crown` is fine for CPU validation but MUST be
     replaced by pre-uploaded resident weights in the GPU loop — never re-upload per
     sub-domain.
  2. **The BaB batch lives on the GPU.** Sub-domain boxes are `gpuArray` columns; the
     whole branch-and-bound loop (split → bound → prune) runs device-side.
  3. **`gather` only the verdicts** (a few booleans / the final answer), never
     intermediate bounds or coefficient matrices.
  4. **Batch with `pagemtimes`/`pagefun`.** The backward CROWN coefficient matrices
     across sub-domains are pages → one batched kernel, not a host-side loop. Use
     `arrayfun` for the elementwise ReLU/MaxPool relaxation (fuses to one CUDA kernel).
  5. **Function overloading is automatic** but only DATA inputs trigger GPU. Keep both
     weights and boxes on the device so every op stays on the GPU.
  6. **Time correctly:** `gputimeit`, or `tic/toc` with `wait(gpuDevice)`; warm up
     once first (JIT/kernel-compile). The MATLAB Profiler cannot time GPU code.

**Validation order:** all soundness validation runs on **CPU** (identical code, host
arrays) — done for IBP and CROWN already. The GPU instance (g4dn, licensed via moving
ENI `eni-0a4349f2e1fa10e5f`) is provisioned only for **throughput measurement** once
the batched form exists.
