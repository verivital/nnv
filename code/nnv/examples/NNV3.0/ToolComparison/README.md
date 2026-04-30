# ToolComparison — NNV vs MathWorks AI Verification Toolbox

Head-to-head comparison of the NNV neural-network verifier against the
MathWorks **AI Verification Toolbox** (formerly the Deep Learning Toolbox
Verification Library). Two complementary halves cover the verification
regimes the paper claims:

- **`acas_rl_tll/`** — refresh of NNV 2.0 / CAV'23 Tables 2–3 on
  feed-forward networks with VNNLIB-style output specifications. NNV's
  star-set reachability is compared against AIVL's
  `estimateNetworkOutputBounds`.
- **`mnist_resnet/`** — first head-to-head on a residual network
  (MNIST-ResNet-8). NNV's star-set reachability is compared against
  AIVL's `verifyNetworkRobustness` with DeepPoly. Enabled by R2024b's
  `additionLayer` support in `verifyNetworkRobustness`.

For each instance we record verdict (verified / violated / unknown /
timeout / error) and per-instance solve time. The schema is uniform
across both halves so a single set of helpers (`tool_utils.m`,
`export_csv.m`, `make_*_table.m`) covers reporting.

## References

- **NNV 3.0 (this work)**: Tumlin, A.M., Manzanas Lopez, D., Johnson,
  T.T., et al. *NNV 3.0 tool paper.* ATVA 2026 (submitted).
- **NNV 2.0 / CAV'23**: Manzanas Lopez, D., Choi, S.W., Tran, H.-D.,
  Johnson, T.T. *NNV 2.0: The neural network verification tool.*
  CAV 2023.
- **MathWorks AI Verification Toolbox**: MathWorks Inc., *AI Verification*,
  R2024a–R2025b. https://www.mathworks.com/help/deeplearning/ai-verification.html

## Algorithm grid (CAV'23-style)

Both halves expose the full NNV 2.0 algorithm grid plus the AIVL APIs
applicable to that half's specification style:

| Tool | Algorithms (acas_rl_tll) | Algorithms (mnist_resnet) |
|------|--------------------------|---------------------------|
| `nnv` | `approx-star`, `relax-star-range-{25,50,75,100}`, `exact-star` | `approx-star`, `relax-star-area-{25,50,75,100}` |
| `mw_estimate` | `estimate-bounds` | (n/a — VNNLIB half-spaces only) |
| `mw_deeppoly` | (n/a) | `deep-poly` |
| `mw_abc` | `alpha-beta-crown` (R2026a bridge required) | (n/a) |

The legacy algorithm string `relax-star-50` is accepted as an alias for
the canonical `relax-star-range-50` (`acas_rl_tll`) or `relax-star-area-50`
(`mnist_resnet`) so existing persisted result files keep matching on
resume.

## Layout

```
NNV3.0/ToolComparison/
├── README.md                          ← this file
├── REPORT.md                          full methodology + paper-ready writeup
├── run_toolcomparison.m               top-level orchestrator (called from run_all.sh)
├── tool_utils.m                       canonical result-row schema + helpers
├── rebuild_for_aivl.m                 strip ScalingLayer / adapter layers for AIVL
├── export_csv.m                       .mat → .csv mirrors for VS Code preview
├── smoke_test.m                       ~30s end-to-end sanity check
│
├── acas_rl_tll/
│   ├── run_acas_rl_tll.m              driver (ACAS Xu / RL / TLLverify)
│   └── results/                       results_<benchmark>.{mat,csv}
│
├── mnist_resnet/
│   ├── run_mnist_resnet.m             driver (MNIST-ResNet-8)
│   ├── train_mnist_resnet.m           native-MATLAB ResNet trainer
│   ├── models/                        trained dlnetwork + 100-image testset
│   └── results/                       expC_<model>.{mat,csv}
│
├── tables/
│   ├── make_acas_rl_tll_table.m       renders Table A
│   ├── make_mnist_resnet_table.m      renders Table C
│   └── out/                           table_A.{tex,txt}, table_C.{tex,txt},
│                                       sanity_report.txt (CAV'23 cross-check)
│
└── scripts/
    └── toolbox_install.m              extracts the AI Verification Toolbox
                                        Support Package tarball into MATLAB userpath
```

## Running

NNV must be installed (`code/nnv/install.m`) and on the MATLAB path. The
runner adds the comparison helpers automatically via `addpath`.

### Default (full grid)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/ToolComparison; run_toolcomparison"
```

Runs both halves with the full CAV'23-style algorithm grid and renders
the paper tables. Wall-clock ~5–6 hours on a 4-core CPU.

### Smoke (~10–15 min)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/ToolComparison; \
    run_toolcomparison('mode','smoke')"
```

5 ACAS networks × {`approx-star`, `relax-star-range-50`} + 5 MNIST
images × {`approx-star`, `relax-star-area-50`}, NNV-only (skips MW
calls). Useful as a CI-style sanity check.

### Half at a time

```matlab
% Only the FC half:
run_acas_rl_tll('benchmarks',{'acas_p3','rl'}, 'algorithms',{'approx-star','exact-star'});

% Only the ResNet half:
run_mnist_resnet('models',{'mnist_resnet8'}, 'epsilons',[1/255]);
```

### Tables only (no verifier compute)

```matlab
make_acas_rl_tll_table;     % from bundled .mat in acas_rl_tll/results/
make_mnist_resnet_table;    % from bundled .mat in mnist_resnet/results/
```

## AI Verification Toolbox setup

The MW-side comparison (`mw_estimate`, `mw_deeppoly`, `mw_abc`) requires
the AI Verification Toolbox Support Package, gated by a MathWorks account
and not redistributable. NNV-only runs work without it; pass
`'tools',{'nnv'}` to skip.

To enable the MW-side, produce a tarball on a workstation that has AIVL
installed. From `~/Documents/MATLAB`:

```sh
tar czf atva26-aivl.tar.gz \
    SupportPackages/R*/toolbox/nnet/supportpackages/aivnv \
    SupportPackages/R*/resources/aivnv \
    SupportPackages/R*/appdata/components/softwaresupportpkgs/nnet/aivnv* \
    SupportPackages/R*/appdata/files/softwaresupportpkgs/nnet/aivnv*
```

Stage the tarball at one of:

- `/home/matlab/addons/atva26-aivl.tar.gz` (NNV3.0 Docker image)
- `/data/aivl/atva26-aivl.tar.gz` (CodeOcean private data mount)
- `<env var TOOLCOMPARISON_AIVL_TARBALL>`
- `scripts/atva26-aivl.tar.gz` (alongside `toolbox_install.m`)

then run:

```matlab
run('scripts/toolbox_install.m');
```

The helper finds the `aivnv` directory wherever it lands (the tarball's
source-release subdir is usually different from the running MATLAB
release; that's fine, AIVL's .m files are version-portable).

## Outputs

Persisted result files are MATLAB tables with the canonical schema
defined in `tool_utils.m`:

```
tool, benchmark, instance_id, status, time, algorithm, timeout, note
```

Each driver writes incrementally and skips instances already present
(`tool_utils.has_instance`), so an interrupted run resumes from where
it left off. CSV mirrors are produced by `export_csv` for VS Code preview.

The bundled `acas_rl_tll/results/` and `mnist_resnet/results/` shipped
with this folder are the **paper-grade results** from the original
ATVA26 standalone artifact (R2025b), produced before the algorithm grid
was extended. New runs that include `approx-star` / `relax-star-*-25` /
etc. add fresh rows alongside without disturbing the existing rows.

## Expected runtime

| Mode | Wall-clock | Notes |
|------|-----------:|-------|
| Smoke | ~10–15 min | 5 ACAS networks + 5 MNIST images, NNV-only |
| acas_rl_tll full (NNV grid only) | ~5 h | Driven by `exact-star` on ACAS Xu (~3–5 h) |
| mnist_resnet full (NNV grid only) | ~2 h | 50 × 4 ε × 5 algorithms |
| Default `run_toolcomparison` | ~6–8 h | Both halves, NNV + AIVL |

The longest pole is `exact-star` on ACAS Xu (CAV'23 reported 10 479 s
for one outlier instance). The default 900 s timeout caps that; raise
it via `'timeout', 14400` to recover all CAV'23 verified cases.

## CAV'23 cross-check

`tables/out/sanity_report.txt` re-renders the comparison against
CAV'23's published exact-star numbers (V=42, X=3, ?=0 for both
ACAS P3 and P4). The cross-check uses (verified + timeout) since
our 900 s cap turns CAV'23's slow-but-verified outliers into timeouts.

## Troubleshooting

**`Inconsistency between affine mapping matrix and dimension of the star
set`** during NNV reach: a TLLverify or RL ONNX file has adapter layers
(VerifyBatchSize, Flatten*) that NNV's matlab2nnv can't strip. The
driver runs the file through `rebuild_for_aivl` first — make sure that
helper is on the path (`addpath` happens automatically when the driver
is invoked, but if you call inner functions directly you'll need to
add it yourself).

**`aivnv:verifyNetwork:DisallowedLayers`** from `verifyNetworkRobustness`:
the network has a `ScalingLayer` or other adapter that AIVL rejects.
The drivers route MW-side calls through `rebuild_for_aivl` for sequential
networks; the residual MNIST-ResNet-8 has no ScalingLayer because it's
trained natively in MATLAB.

**`verifyNetworkRobustness` fails on additionLayer**: requires R2024b+.
On older releases, the MW-side rows for `mnist_resnet` will record
`error`; the bundled paper-grade results in `mnist_resnet/results/`
are still rendered into Table C by `make_mnist_resnet_table`.
