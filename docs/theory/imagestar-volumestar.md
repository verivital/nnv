# ImageStar & VolumeStar

## ImageStar

The **ImageStar** extends the Star set concept to multi-channel 2D images, preserving spatial structure through convolutional layers.

### Definition

$$\Theta = \{ x \in \mathbb{R}^{H \times W \times C} \mid x = c + \sum_{i=1}^{m} \alpha_i v_i, \;\; C\alpha \leq d \}$$

where:
- $c \in \mathbb{R}^{H \times W \times C}$ is the center image
- $v_i \in \mathbb{R}^{H \times W \times C}$ are generator images
- $\alpha \in \mathbb{R}^m$ are predicate variables
- $C\alpha \leq d$ are linear constraints

### Key Insight

The ImageStar preserves the **spatial structure** of images, which is critical for efficient convolution. Rather than flattening to a 1D Star (which would make convolution a dense matrix multiply), the ImageStar applies convolution directly to the center and each generator image:

$$\text{Conv}(\Theta) = \{ y \mid y = \text{Conv}(c) + \sum_i \alpha_i \cdot \text{Conv}(v_i), \;\; C\alpha \leq d \}$$

This means the convolution of an ImageStar is computed as:
1. Convolve the center image
2. Convolve each generator image independently
3. The constraints remain unchanged

This is vastly more efficient than converting to a dense Star representation.

### Operations

| Operation | ImageStar Behavior |
|-----------|-------------------|
| Convolution (Conv2D) | Apply conv to center and each generator |
| Batch normalization | Scale and shift center and generators |
| Average pooling | Apply pooling to center and generators |
| Max pooling | Requires LP-based splitting or approximation |
| ReLU | Per-pixel exact or approximate analysis |
| Flatten | Convert to 1D Star set |

## VolumeStar

The **VolumeStar** further extends ImageStar to 4-dimensional data: 3D spatial volumes or video frames.

### Definition

$$V = \{ x \in \mathbb{R}^{H \times W \times C \times F} \mid x = c + \sum_{i=1}^{m} \alpha_i v_i, \;\; P\alpha \leq q \}$$

where:
- $c \in \mathbb{R}^{H \times W \times C \times F}$ is the center volume
- $v_i \in \mathbb{R}^{H \times W \times C \times F}$ are generator volumes
- The 4th dimension $F$ can represent depth (for 3D medical images) or frames (for video)

### Supported Layers

The VolumeStar propagates through:
- **Conv3D**: 3D convolution applied to center and generators
- **AveragePooling3D**: 3D spatial pooling
- **ReLU / activation layers**: Per-voxel analysis
- **Flatten**: Convert to 1D Star

### Applications

- **Video classification**: Verify C3D/I3D networks under temporal and spatial perturbations
- **3D medical imaging**: Verify classifiers on MRI/CT volumes under noise and artifact perturbations
- **Spatio-temporal robustness**: Perturbations that affect both space and time dimensions

### Hierarchy

The set representations form a natural hierarchy:

```
Star (1D vectors)
 └── ImageStar (2D: H × W × C)
      └── VolumeStar (3D: H × W × C × F)

Zono (1D zonotopes)
 └── ImageZono (2D zonotopes)

GraphStar (graph-structured: N × F)
```

Each level specializes the parent representation for higher-dimensional data while maintaining the same constraint structure ($C\alpha \leq d$).

## References

- H.-D. Tran, S. Bak, W. Xiang, T.T. Johnson, "Towards Verification of Large Convolutional Neural Networks Using ImageStars," CAV 2020
- S. Sasaki, D. Manzanas Lopez, P.K. Robinette, T.T. Johnson, "Robustness Verification of Video Classification Neural Networks," FormaliSE 2025
