# ToolComparison — NNV vs MathWorks AI Verification Library

Head-to-head comparison of NNV against the MathWorks **AI Verification Library**
(AIVL; formerly Deep Learning Toolbox Verification Library). Two halves cover
the verification regimes the NNV3.0 paper claims:

- **`acas_rl_tll/`** — feed-forward and CNN networks with VNNLIB-style output
  specifications. NNV's star-set reachability is compared against AIVL's
  `estimateNetworkOutputBounds`. Active benchmarks:
  - ACAS Xu p3, p4 (CAV'23 baseline; full algorithm grid)
  - RL controllers (cartpole, lunarlander, ...; CAV'23, `approx-star`)
  - **OVAL21** — 3 CIFAR CNNs × 10 properties (VNNCOMP 2021)
  - **Collins RUL** — 62 properties across 3 RUL CNNs (VNNCOMP 2022)

  TLLverify (CAV'23) is retired; bundled rows preserved at
  [`acas_rl_tll/legacy/`](acas_rl_tll/legacy/).
- **`mnist_resnet/`** — first head-to-head on a residual network
  (MNIST-ResNet-8). NNV's star-set reachability vs AIVL's
  `verifyNetworkRobustness` with DeepPoly. Enabled by R2024b's
  `additionLayer` support.

For each instance we record verdict (verified / violated / unknown / timeout
/ error) and per-instance solve time. Schema is uniform across both halves:
a MATLAB `table` with columns
`{tool, benchmark, instance_id, status, time, algorithm, timeout, note}`.

## References

- **NNV 3.0 (this work)**: Tumlin, A.M., Manzanas Lopez, D., Johnson, T.T.,
  et al. *NNV 3.0 tool paper.* ATVA 2026 (submitted).
- **NNV 2.0 / CAV'23**: Manzanas Lopez, D., Choi, S.W., Tran, H.-D.,
  Johnson, T.T. *NNV 2.0: The neural network verification tool.* CAV 2023.
- **MathWorks AI Verification Library**: MathWorks Inc., *AI Verification*,
  R2024a–R2025b. https://www.mathworks.com/help/deeplearning/ai-verification.html

## Algorithm grid (per-benchmark, NNV 2.0 / CAV'23-aligned)

Each benchmark uses the algorithm grid from its original paper. ACAS Xu
p3/p4 use the full CAV'23 grid; RL is `approx-star` only per NNV 2.0
methodology; the new VNNCOMP-derived FC/CNN benchmarks extend this story.
For MNIST-ResNet-8 we use a single `relax-star-area-50` configuration
that matches the residual-network methodology.

| Tool | acas_p3 / p4 | rl | oval21 | collins_rul | mnist_resnet8 |
|------|--------------|----|--------|-------------|---------------|
| `nnv` | approx-star + relax-star-range-{25,50,75,100} + exact-star | approx-star | approx-star + relax-star-area-{25,50,75,100} | approx-star + relax-star-area-{25,50,75,100} + exact-star | relax-star-area-50 |
| `mw_estimate` | estimate-bounds | estimate-bounds | estimate-bounds | estimate-bounds | (n/a — argmax robustness) |
| `mw_deeppoly` | (n/a — VNNLIB half-spaces only) | (n/a) | (n/a) | (n/a) | deep-poly |
| `mw_abc` | alpha-beta-crown (R2026a bridge required) | (n/a) | (n/a) | (n/a) | (n/a) |

The driver `algorithms_for(tool, benchmark)` enforces these defaults; pass
`'algorithms', {...}` to override. The legacy algorithm string `relax-star-50`
is accepted as an alias for `relax-star-range-50` (FC) or `relax-star-area-50`
(CNN), so persisted result rows from older runs keep matching on resume.

## Headline results

Sanity check: `tables/out/sanity_report.txt` confirms ACAS p3/p4 match
CAV'23 exact-star (`Verified + Timeout = 42` exactly).

| Benchmark | NNV best | AIVL | Highlight |
|-----------|----------|------|-----------|
| acas_p3 | exact-star V=32, X=3, T=10 (129 s) | estimate-bounds V=2, X=3, ?=40 (0.05 s) | NNV verifies 16× more |
| acas_p4 | exact-star V=39, X=3, T=3 (170 s) | estimate-bounds V=0, X=1, ?=44 (0.03 s) | NNV verifies 39 of 45 |
| rl | approx-star V=32, X=14, ?=4 (0.09 s) | estimate-bounds V=20, X=11, ?=19 (0.02 s) | NNV verifies 60% more |
| oval21 | approx-star X=10 (5 timeouts, 102 s) | estimate-bounds X=9 (0.18 s) | NNV finds 1 more counterexample |
| collins_rul | approx-star V=10, X=47, ?=5 (2.3 s) | estimate-bounds V=10, X=47, ?=5 (0.17 s) | identical verdicts |
| mnist_resnet8 | relax-star-area-50 V=50/50 (6–39 s, ε=1/255–8/255) | deep-poly V=50/50 (~43 s flat) | NNV 3.1×–7.4× faster at small ε |

See the ATVA 2026 paper for full per-algorithm tables.

## Layout

```
NNV3.0/ToolComparison/
├── README.md                          ← this file
├── run_toolcomparison.m               top-level orchestrator
│
├── utils/                             helpers shared across the suite
│   ├── tool_utils.m                   canonical result-row schema + helpers
│   ├── rebuild_for_aivl.m             strip ScalingLayer / adapter layers for AIVL
│   ├── export_csv.m                   .mat → .csv mirrors for VS Code preview
│   ├── smoke_test.m                   ~30 s end-to-end sanity check
│   ├── toolbox_install.m              extract AIVL Support Package tarball
│   └── canonicalize_bundled_results.m  one-shot rename of legacy `relax-star-50` rows
│
├── acas_rl_tll/                       FC + image-VNNLIB half
│   ├── run_acas_rl_tll.m              driver (ACAS, RL, oval21, collins_rul)
│   ├── collins_rul_cnn_2022/          VNNCOMP'22 RUL CNN assets (3 ONNX + 41 vnnlib)
│   ├── results/                       results_<benchmark>.{mat,csv}
│   └── legacy/                        retired benchmarks (TLLverify; see below)
│
├── mnist_resnet/                      argmax-robustness half
│   ├── run_mnist_resnet.m             driver
│   ├── train_mnist_resnet.m           native-MATLAB ResNet trainer
│   ├── models/                        trained dlnetworks + testsets
│   └── results/                       expC_<model>.{mat,csv}
│
└── tables/
    ├── make_acas_rl_tll_table.m       renders paper Tables 5 (FC) + 6 (CNN)
    ├── make_mnist_resnet_table.m      renders paper Table 7 (MNIST-ResNet-8)
    └── out/                           {table_A,table_C,sanity_report}.{tex,txt}
```

OVAL21 and ACAS / RL ONNX assets live in
[`code/nnv/examples/NNV2.0/Submission/CAV2023/NNV_vs_MATLAB/`](../../NNV2.0/Submission/CAV2023/NNV_vs_MATLAB/),
referenced by the driver (no copy needed).

## Running

NNV must be installed (`code/nnv/install.m`) and on the MATLAB path. The
runner adds the comparison helpers automatically via `addpath`.

### Default (full grid)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/ToolComparison; run_toolcomparison"
```

Runs both halves with the per-benchmark algorithm grids. Wall-clock ~3–5 h
on a 4-core CPU, dominated by `exact-star` on ACAS Xu (≤900 s/instance).

### Smoke (~10–15 min)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/ToolComparison; \
    run_toolcomparison('mode','smoke')"
```

### Half at a time

```matlab
% FC half:
run_acas_rl_tll('benchmarks',{'acas_p3','rl'}, 'algorithms',{'approx-star','exact-star'});
% ResNet half:
run_mnist_resnet('models',{'mnist_resnet8'}, 'epsilons',[1/255]);
```

### Tables only (no verifier compute)

```matlab
make_acas_rl_tll_table;     % from bundled .mat in acas_rl_tll/results/
make_mnist_resnet_table;    % from bundled .mat in mnist_resnet/results/
```

## AI Verification Library setup

The MW-side comparison (`mw_estimate`, `mw_deeppoly`, `mw_abc`) requires
the AIVL Support Package, gated by a MathWorks account and not
redistributable. NNV-only runs work without it; pass `'tools',{'nnv'}`
to skip.

To enable the MW-side, produce a tarball on a workstation that has AIVL
installed. From `~/Documents/MATLAB`:

```sh
tar czf atva26-aivl.tar.gz \
    SupportPackages/R*/toolbox/nnet/supportpackages/aivnv \
    SupportPackages/R*/resources/aivnv \
    SupportPackages/R*/appdata/components/softwaresupportpkgs/nnet/aivnv* \
    SupportPackages/R*/appdata/files/softwaresupportpkgs/nnet/aivnv*
```

Stage at one of `/home/matlab/addons/atva26-aivl.tar.gz` (Docker),
the env var `TOOLCOMPARISON_AIVL_TARBALL`, or `utils/atva26-aivl.tar.gz`.
Then run:

```matlab
run('utils/toolbox_install.m');
```

The helper finds the `aivnv` directory wherever it lands, so a tarball
made on R2025b extracts cleanly into R2024a etc. AIVL is not officially
registered as an installed Support Package; the `addpath` in
`startup.m` is what makes `verifyNetworkRobustness` /
`estimateNetworkOutputBounds` resolve to the bundled implementation
rather than the matlabroot stub. Worker-side `addpath` is repeated in
the verifier's MW path because parfeval workers do not auto-run
`startup.m`.

## Infrastructure notes

**ScalingLayer dispatch (R2025a+ ONNX importer).** R2025a wraps each FC's
bias-add as `nnet.cnn.layer.ScalingLayer` instead of
`nnet.onnx.layer.ElementwiseAffineLayer`. Two patches in NNV's engine
(both upstream-ready):

- [`code/nnv/engine/nn/layers/ElementwiseAffineLayer.m`](../../../engine/nn/layers/ElementwiseAffineLayer.m) —
  `parse()` accepts both class names.
- [`code/nnv/engine/utils/matlab2nnv.m`](../../../engine/utils/matlab2nnv.m) —
  dispatch updated.

**`rebuild_for_aivl`** ([utils/rebuild_for_aivl.m](utils/rebuild_for_aivl.m)).
AIVL rejects `ScalingLayer` and `VerifyBatchSizeLayer`. The helper folds
ScalingLayer into the preceding FC's bias for FC nets and surgically
strips disallowed adapter layers for CNN nets.

**GLPK MEX.** NNV's [engine/utils/lpsolver.m](../../../engine/utils/lpsolver.m)
falls back to GLPK on degenerate LPs. The Docker `tbxmanager install
glpkmex` step had been silently failing in the build; the bind-mounted
toolbox path makes GLPK visible to the container at runtime.

**VNNLIB output spec parser.** NNV's `load_vnnlib_matlab.m` has a
paren-counting bug on `(or (and (<= Y_i Y_j)))` outputs (oval21,
collins_rul). The MW-side `verifyImageVnnlib` bypasses the parser by
loading the property via `load_vnnlib` (works) then evaluating the
disjunctive halfspaces against AIVL output bounds via `mw_vnnlib_status`.

## Outputs

Persisted result files are MATLAB tables with the canonical schema in
[`utils/tool_utils.m`](utils/tool_utils.m):

```
tool, benchmark, instance_id, status, time, algorithm, timeout, note
```

Drivers write incrementally and skip instances already present
(`tool_utils.has_instance`), so an interrupted run resumes from where it
left off. CSV mirrors are produced by `export_csv` for VS Code preview.

The bundled `acas_rl_tll/results/` and `mnist_resnet/results/` shipped
with this folder hold the full Tables 5, 6, and 7 results behind the paper.
For ACAS p3/p4 and RL, paper-grade rows carried over from the CAV'23
NNV 2.0 artifact (re-run on R2025b) are preserved; subsequent runs from
this tree added the extended algorithm grid alongside without disturbing
them. OVAL21 and Collins RUL are fully fresh from this tree.

## Expected runtime

| Mode | Wall-clock | Notes |
|------|-----------:|-------|
| Smoke | ~10–15 min | 5 ACAS networks + 5 MNIST images, NNV-only |
| FC half (NNV grid only) | ~3–5 h | ACAS Xu full grid drives the runtime |
| ResNet half (NNV only) | ~10 min | 50 × 4 ε × 1 algorithm (`relax-star-area-50`) |
| Default `run_toolcomparison` | ~3–5 h | Both halves, NNV + AIVL |

The longest pole is `exact-star` on ACAS Xu (CAV'23 reported 10 479 s
for one outlier). The default 900 s timeout caps that; raise it via
`'timeout', 14400` to recover all CAV'23 verified cases. OVAL21's
`cifar_wide_kw` at large ε is the second-longest pole (~250–500 s per
instance with `relax-star-area-25`).

