# NNV 3.0: Major Contributions and New Features

## Overview

NNV 3.0 represents a significant evolution of the Neural Network Verification tool, building on the foundations laid in NNV 1.0 (CAV 2020) and NNV 2.0 (CAV 2023). While previous versions focused on expanding architectural support, NNV 3.0 introduces methodologies for verifying new classes of properties and handling previously intractable data formats.

This document summarizes the major contributions in NNV 3.0 compared to earlier versions.

---

## Version History

### NNV 1.0 (CAV 2020)
- Initial release with Star-based reachability analysis
- Support for fully-connected and convolutional neural networks (FFNNs, CNNs)
- ReLU activation function verification
- ImageStar representation for image inputs
- Neural Network Control Systems (NNCS) verification

### NNV 2.0 (CAV 2023)
- Recurrent Neural Networks (RNNs) support
- Semantic Segmentation Neural Networks (SSNNs)
- Neural Ordinary Differential Equations (Neural ODEs)
- ONNX and VNNLIB format support
- VNN-COMP participation
- Relaxed reachability methods for scalability

### NNV 3.0 (Current Release)
Major new features detailed below.

---

## NNV 3.0 Major Research Contributions

### 1. VolumeStar: Spatio-Temporal and Volumetric Data Verification

**First-of-its-kind verification for video and 3D data.**

A major challenge in neural network verification is handling high-dimensional data, particularly data with both spatial and temporal components like videos. NNV 3.0 introduces the VolumeStar, a new set representation specifically designed for verifying models that operate on 3-dimensional data.

- Novel set representation extending Star/ImageStar to 4D (H × W × C × F)
- Enables robustness verification of video classification networks (C3D, I3D)
- Supports volumetric medical images (3D MRI, CT scans)
- Handles spatial AND temporal perturbations
- Compact representation for efficient layer-by-layer analysis

The VolumeStar is defined as:

$$V = \{ x \in \mathbb{R}^{H \times W \times C \times F} \mid x = c + \sum_{i=1}^{m} \alpha_i v_i, \text{s.t. } P\alpha \le q \}$$

where c is the center "volume" and v_i are the generator "volumes" representing the basis of variation.

### 2. FairNNV: Fairness Verification Module

**Formal verification of algorithmic fairness.**

As AI models are increasingly deployed for high-stakes decisions (e.g., loan approvals, credit scoring), ensuring their fairness has become a critical requirement. FairNNV is a module in NNV 3.0 that enables formal verification of algorithmic fairness via reachability analysis.

- **Counterfactual Fairness**: Proves predictions unchanged when sensitive attributes (gender, race) are altered
- **Individual Fairness**: Verifies similar individuals receive similar predictions
- **Verified Fairness (VF) Score**: Quantitative, formally grounded measure
- Reachability-based approach guarantees hold for ALL inputs, not just samples

Unlike statistical auditing which considers only sampled inputs, FairNNV soundly verifies fairness over continuous input regions, providing guarantees that hold for all inputs within specified bounds.

### 3. Probabilistic Verification via Conformal Prediction

**Model-agnostic, scalable verification for large networks.**

Some problems such as semantic segmentation are inherently expensive for neural network verification as the output reachable set being computed is intractable given scalability of existing techniques. NNV 3.0 includes a response to these challenges via a data-driven probabilistic reachability algorithm based on conformal inference.

- Data-driven reachability based on conformal inference
- Works as drop-in alternative when full worst-case guarantees are infeasible
- User specifies: perturbation space, sampling distribution, miscoverage level
- Returns probabilistic output sets with formal coverage guarantees
- Particularly effective for semantic segmentation and large architectures

### 4. ModelStar: Weight Perturbation Analysis

**Verification under parameter uncertainties (quantization, hardware errors).**

Parameter perturbations pose safety-critical concerns for model deployment after training. These perturbations may be caused by device variations as well as quantization. NNV 3.0 extends the framework with star set augmentations that capture interval-based parameter perturbations for linear layers.

- Star set augmentations for interval-based parameter perturbations
- **Single-layer perturbation**: Exact star set representation
- **Multi-layer perturbation**: Sound over-approximation
- Implemented for fully-connected and 2D convolutional layers
- Addresses deployment safety: quantization uncertainty, device variations

### 5. Time-Dependent Neural Networks

**Variable-length time series verification.**

While NNV 2.0 introduced support for RNNs, the analysis was primarily focused on fixed-length input sequences. NNV 3.0 incorporates a new methodology for verifying time-dependent neural networks that can process variable-length sequences.

- Extends RNN verification beyond fixed-length sequences
- Handles time horizon T ∈ [T_min, T_max]
- Computes union of reachable sets over possible time horizons
- Applications: NLP, sensor data analysis, sequential decision systems

### 6. Malware Detection Benchmark

**New verification domain: cybersecurity.**

The use of DNNs for malware detection has shown great promise, but these models are also vulnerable to adversarial attacks where malware authors make small modifications to evade detection. NNV 3.0 introduces a new benchmark for verifying the robustness of such models.

- Robustness verification of DNN-based malware classifiers
- Feature-space perturbations mimicking adversarial evasion
- Formal proof that modifications cannot cause misclassification
- New benchmark for community evaluation

---

## Software Engineering Improvements

### Enhanced Usability

- **NNVVERSION() function**: Programmatic version querying
- **MATLAB R2023a+ version check**: Clear compatibility requirements
- **Three-tier dependency checking**: Critical/Important/Optional classification
- **check_nnv_setup.m**: Quick diagnostic tool for troubleshooting
- **"Ready!" confirmation**: Clear startup success indicator

### CI/CD Pipeline

- GitHub Actions for automated testing
- Every commit triggers unit tests and regression tests
- Ensures contributions don't break existing functionality
- Stable platform for community contributions

### Documentation and Tutorials

- **QuickStart folder**: test_installation.m, simple_verification.m
- **Tutorial series**: ACAS Xu, Image Classifiers, NNCS
- **MATLAB Online support**: Try NNV without local installation
- Conference tutorials: EMSOFT, SPIE, DSN, IAVVC

### Test Infrastructure

- 833 tests (full suite), 470 tests (quick mode) - 100% pass rate
- 48.6% code coverage (86.2% nn/, 88.9% set/)
- Platform-aware installation (Windows, Linux, macOS)
- Regression testing with baseline comparison

### Bug Fixes

- lpsolver.m version checking compatibility with newer MATLAB
- BatchNormalizationLayer constructor (name-value pairs)
- install.m typos and messaging improvements

---

## Comparison: NNV vs Other Verification Tools

| Architecture/Application | NNV 3.0 | α,β-Crown | CORA | Marabou | nnenum |
|-------------------------|---------|-----------|------|---------|--------|
| FFNN | ✅ | ✅ | ✅ | ✅ | ✅ |
| CNN | ✅ | ✅ | ✅ | ✅ | ✅ |
| RNN | ✅ | ✅ | ❌ | ❌ | ❌ |
| SSNN | ✅ | ~ | ❌ | ❌ | ❌ |
| Time-dep. NNs | ✅ | ❌ | ❌ | ❌ | ❌ |
| Video NNs | ✅ | ❌ | ❌ | ❌ | ❌ |
| 3D (volumetric) | ✅ | ❌ | ❌ | ❌ | ❌ |
| Neural ODEs | ✅ | ❌ | ❌ | ❌ | ❌ |
| NNCS | ✅ | ❌ | ✅ | ❌ | ❌ |
| Fairness | ✅ | ❌ | ❌ | ❌ | ❌ |
| Weight Perturbation | ✅ | ❌ | ❌ | ❌ | ❌ |

NNV 3.0 is the most comprehensive verification framework supporting the widest range of architectures and verification properties.

---

## Getting Started with NNV 3.0

### 1. Verify Installation

```matlab
check_nnv_setup()
```

### 2. Check Version

```matlab
NNVVERSION()  % Returns 'NNV v3.0.0'
```

### 3. Run QuickStart

```matlab
cd examples/QuickStart
test_installation
simple_verification
```

### 4. Verify ONNX + VNNLIB

```matlab
% Load network from ONNX
dlnetwork = importNetworkFromONNX(path_to_net, InputDataFormats='BCSS');
net = matlab2nnv(dlnetwork);

% Define reachability options
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

% Verify property
res = net.verify_vnnlib(vnnlib_file, reachOptions);
```

---

## Citation

If you use NNV in your research, please cite:

### NNV 1.0 (CAV 2020)
```bibtex
@inproceedings{tran2020nnv,
  author = {Tran, Hoang-Dung and Yang, Xiaodong and Manzanas Lopez, Diego and
            Musau, Patrick and Nguyen, Luan Viet and Xiang, Weiming and
            Bak, Stanley and Johnson, Taylor T.},
  title = {NNV: The Neural Network Verification Tool for Deep Neural Networks
           and Learning-Enabled Cyber-Physical Systems},
  booktitle = {Computer Aided Verification (CAV)},
  year = {2020},
  doi = {10.1007/978-3-030-53288-8_1}
}
```

### NNV 2.0 (CAV 2023)
```bibtex
@inproceedings{lopez2023nnv,
  author = {Lopez, Diego Manzanas and Choi, Sung Woo and
            Tran, Hoang-Dung and Johnson, Taylor T.},
  title = {NNV 2.0: The Neural Network Verification Tool},
  booktitle = {Computer Aided Verification (CAV)},
  year = {2023},
  doi = {10.1007/978-3-031-37703-7_19}
}
```

### NNV 3.0 (FM 2026)
*Paper in preparation*

---

## Related Publications

- **FairNNV**: Tumlin et al., "FairNNV: The Neural Network Verification Tool For Certifying Fairness", ICAIF 2024
- **Time-Dependent Networks**: Pal et al., "Robustness Verification of Deep Neural Networks using Star-Based Reachability Analysis with Variable-Length Time Series Input", FMICS 2023
- **Malware Detection**: Robinette et al., "Case Study: Neural Network Malware Detection Verification", FormaliSE 2024
- **Semantic Segmentation**: Tran et al., "Robustness Verification of Semantic Segmentation Neural Networks using Relaxed Reachability", CAV 2021

---

## Contact

For questions or issues, please visit: [github.com/verivital/nnv/issues](https://github.com/verivital/nnv/issues)
