# Sound star-set verification of the VNN-COMP 2023 ViT benchmark

This example verifies the **VNN-COMP 2023 Vision Transformer** track
(`shizhouxing/ViT_vnncomp2023`) entirely on NNV `Star` sets, using the sound
softmax-attention reachability primitives added in
[`engine/nn/funcs/SoftmaxAttn.m`](../../../engine/nn/funcs/SoftmaxAttn.m).

- **Models** (CIFAR-10, 3×32×32): `pgd_2_3_16` (patch 16, depth 2, 5 tokens) and
  `ibp_3_3_8` (patch 8, depth 3, 17 tokens). Both: embed dim 48, 3 heads, head
  dim 16, MLP hidden 96, scale 1/4.
- **Property**: L∞ robustness, ε = 1/255 (normalised), argmax preservation.
- **Architecture** (`ViT_BN`, mirrors the exported ONNX): Conv patch-embed →
  `[cls]` token → `+pos` → *depth* × ( `x = x + Attn(BN(x))`; `x = x + FF_ReLU(BN(x))` )
  → ReduceMean over tokens → BN → Linear head. The export replaced training
  LayerNorm with **BatchNorm** (a fixed affine map in eval) and pools with
  ReduceMean, so the **only nonlinearities are softmax attention and the FF ReLU**.

## The sound attention reach (`SoftmaxAttn`)

The benchmark's attention is `softmax(scale · Q·Kᵀ)·V` with uncertainty on Q, K,
*and* V — a bilinear-then-softmax-then-bilinear composition that NNV could not
previously bound soundly (`DynamicMatmulLayer.reach` errored;
`ScaledDotProductAttentionLayer` box-lifted and failed loud on multi-token input).
`SoftmaxAttn` provides the missing pieces, each Monte-Carlo–validated against the
`linprog` containment oracle (see `tests/nn/attention/`):

| method | what it bounds | soundness |
|---|---|---|
| `intervalMatMul` / `bilinearMatMulStar` | set@set `A·B` (Q·Kᵀ, A·V) | Rump midpoint–radius interval matmul; encloses every product |
| `correlatedRowSoftmaxBounds` | `softmax(S)` over a logit box | **exact** per-element range via the monotonicity corners |
| `avEnvelopeStar` | symbolic `A·V`, A interval, V a Star | sign-aware McCormick; affine in V's predicates (value path stays correlated) |
| `prefixAdd` / `prefixConcat` | residual add / head assembly | provenance-safe prefix alignment (no structural-equality under-approx) |
| `singleHeadAttn` / `selfAttentionReach` | one attention sublayer | composition of the above |

These are wired into `DynamicMatmulLayer` (sound interval matmul when operand
shapes are set) and `ScaledDotProductAttentionLayer` (sound multi-token attention
across all star-family methods), so the toolbox layer ecosystem gains sound
attention too — not just this driver.

## Running it

```matlab
% 1. one-time: generate the weight/instance bundles from the benchmark ONNX
%    (needs python with onnx, onnxruntime, scipy; see extract_weights.py header)
%    > python extract_weights.py        % writes pgd_2_3_16.mat, ibp_3_3_8.mat

% 2. in MATLAB (with engine on the path):
addpath(genpath('../../../engine'));
M = ViTReach.load('ibp_3_3_8.mat');
img = squeeze(M.images(1,:,:,:));               % raw [0,1] image
[lb, ub] = ViTReach.epsBox(M, img, 1/255);      % normalised eps-box
[robust, margins] = ViTReach.verify(M, lb, ub, M.labels(1));   % 1 robust, 2 unknown

% sweep + certified radius:
run_benchmark('ibp_3_3_8', 15, true);

% sound branch-and-bound (input split):
[status, info] = ViTReach.verifyBaB(M, img, M.labels(1), 1/255, struct('babMaxNodes',50));
```

`ViTReach.evaluate` is a faithful forward pass (parity to onnxruntime ≈ 1.8e-6,
locked by `test_ViTReach_parity.m`), so the reach verifies the *actual* deployed
ONNX function.

## Results and the role of branch-and-bound

Single-shot sound reach is fast (≈0.8 s pgd, ≈2.9 s ibp) and **sound** (0 sample
escapes across the test suite). At the full ε = 1/255 it certifies essentially
**0/200** — *as expected and by construction*: the benchmark's
`generate_properties.py` keeps only images where PGD fails **and** vanilla CROWN
cannot certify, i.e. exactly the instances every *incomplete* (non-BaB) method
fails on. The published sound score (α,β-CROWN **79/200**) comes entirely from
branch-and-bound.

What the sound reach *does* show: the `ibp_3_3_8` margins approach 0 at full ε
(≈ −0.5, not vacuous) and it **certifies a real radius ε\* ≈ 0.45–0.56 /255**
(median), i.e. sound robustness at roughly half the benchmark ε. The `pgd_2_3_16`
model is hostile to interval bounds (margins ≈ −1000, ε\* ≈ 0), matching
auto_LiRPA IBP.

`ViTReach.verifyBaB` adds a **sound input-splitting BaB** driven by these layers.
Input BaB is theoretically weak in 3072-d (one-pixel splits barely tighten — which
is *why* α,β-CROWN uses neuron/ReLU-split β-CROWN), so it does not close the gap to
1/255; it is included as a sound, working BaB integrated with the new layers and a
template for the neuron-split extension. See `STATUS.md` for measured numbers and
the β-CROWN integration path.
