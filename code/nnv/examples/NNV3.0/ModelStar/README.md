# ModelStar — Weight-Perturbation Verification on an MNIST MLP

Sound reachability-based verification under **weight perturbations**:
for each fully-connected layer of an MNIST classifier, the verifier
sweeps a perturbation magnitude *f* (as a fraction of that layer's
weight-range) and asks: across *all* perturbations of that magnitude,
does the network still classify each held-out image correctly?

The result is a per-layer × per-fraction count of how many of 100 test
images stay safely classified. The technique extends NNV's Star sets
to capture interval-bounded weight perturbations, computing **exact**
reachable sets for single-layer perturbations and sound
over-approximations for multi-layer cases.

## References

- **ModelStar (CAI 2025)**: Zubair, M.U., Johnson, T.T., Basu, K.,
  Abbas, W. *Verification of neural network robustness against weight
  perturbations using star sets.* IEEE Conference on Artificial
  Intelligence (CAI) 2025.
- **ModelStar (extended, JAIR 2026)**: Zubair, M.U., Johnson, T.T.,
  Basu, K., Abbas, W. *ModelStar: Reachability analysis-based safety
  verification of neural networks against model perturbations.* JAIR
  special issue *Integration of Logical Constraints in Deep Learning*
  (2026, accepted).

## Layout

```
NNV3.0/ModelStar/
├── README.md
├── run_expt_for_compute.m   Top-level runner (entry point)
├── conv_expt_any_layer.m    Reachability + verification per layer × fraction
├── EXPT.m                   Class managing experiment config + plotting
├── build_template.m         Programmatic empty-template builder (replaces YAML)
├── mnist_model_fc.mat       MNIST MLP weights (fc_1 .. fc_6)
├── results/                 Output: MNIST_MLP.mat (populated config + verdicts)
└── runtime/                 MATLAB diary log
```

The MNIST MLP is the same 5-hidden-layer network the paper evaluates
(layer sizes 1024, 512, 256, 256, 256). The default sweep verifies the
last three FC layers — `fc_4`, `fc_5`, `fc_6` — against the same
perturbation grid as paper Table 3.

## Running

NNV must be installed (`code/nnv/install.m`) and on the MATLAB path
inside the Docker image.

### Default sweep (paper-aligned)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/ModelStar; run_expt_for_compute"
```

Sweeps fc_4, fc_5, fc_6 × 4 fractions each × 100 images, comparing
ModelStar against the Certificated-Robust ("Towards") and
Formal-Robust ("Formalizing") baselines.

### Smoke (single layer)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/ModelStar; \
    n_layers_to_run_for_from_yaml_file = 1; run_expt_for_compute"
```

Runs only `fc_6` (4 fractions × 100 images).

## Configuration parameters

Edit the layer × fraction grid in
[`build_template.m`](build_template.m). The currently shipped grid
matches the paper:

| Layer | Fractions (% of layer weight range) |
|-------|--------------------------------------|
| fc_6  | 0.5, 1.0, 1.5, 2.0  |
| fc_5  | 0.1, 0.2, 0.3, 0.4  |
| fc_4  | 0.1, 0.2, 0.3, 0.4  |
| fc_3  | 0.05, 0.1, 0.15, 0.2 (declared but not in default sweep) |
| fc_2  | 0.04, 0.06, 0.08, 0.1 (declared but not in default sweep) |
| fc_1  | 0.02, 0.03, 0.04, 0.05 (declared but not in default sweep) |

Top-of-runner overrides:

| Variable | Default | Notes |
|---|---|---|
| `n_layers_to_run_for_from_yaml_file` | `3` | Number of layers from the *end* (fc_6, fc_5, fc_4 with default = 3) |
| `data.n_images`                      | `100` | Set inside `build_template.m` |

## Outputs

- `results/MNIST_MLP.mat` — populated experiment struct. Each
  `data.layers{k}.fracs{m}.percents` stores the {ModelStar,
  Cert-Robust, Formal-Robust} verified counts for that layer × fraction
  combination, and `.times` stores per-sample timing.
- `runtime/conv_expt_any_layer.log` — MATLAB diary.

To view a summary plot post-run:

```matlab
e = EXPT('ModelStar/results/MNIST_MLP'); e.plot_results
```

This folder generates **no figures by default** — figures are produced
on demand from `EXPT.plot_results`.

## Expected runtime

Measured on an NVIDIA RTX 4060 host running the `nnv3.0` Docker image
(MATLAB R2024b), CPU verification:

- **Smoke** (1 layer = fc_6, 4 fractions × 100 images): **~1 minute**
- **Default sweep** (3 layers × 4 fractions × 100 images): **~70 minutes**

The `fc_4` layer dominates the runtime (~30+ minutes) because deeper
perturbed layers compose with more downstream affine + ReLU
operations.
