# OVAL21 Benchmarks.

This benchmark set includes 3 ReLU-based convolutional networks, provided in the `nets` folder,
 which were robustly trained using [1] against l_inf perturbations of radius eps=2/255 on CIFAR10.

The verification properties, provided in the `vnnlib` folder, represent adversarial
robustness to infinity norm perturbations whose radius has been determined individually
for each image.

The same networks and similary-derived perturbation radii were employed in a number of recent works: 
[A](https://openreview.net/forum?id=B1evfa4tPB), [B](http://proceedings.mlr.press/v124/bunel20a/bunel20a.pdf), [C](https://openreview.net/forum?id=uQfOy7LrlTR), [D](https://openreview.net/forum?id=nVZtXBI6LNn), [E](https://arxiv.org/abs/2103.06624), [F](https://arxiv.org/abs/2105.14644).

## Dataset details

Here, we outline the procedure employed to create the dataset. *In short, we rely on commonly employed lower and upper 
bounds to the adversarial loss to exclude perturbation radii that yield trivial properties*. 

10 correctly classified images are randomly sampled from the entire CIFAR10 test set.
For each image, we look for a perturbation radius (epsilon) in the range [0, 16/255].
First, a binary search is run to find the largest epsilon value for which a popular iterative adversarial attack [2] 
fails to find an adversarial example. Then, a second binary search is run to find the smallest epsilon value
for which bounds from the element-wise convex hull of the activations (with fixed intermediate bounds from [1, 4]) [3] fail to prove robustness.
Both binary search procedures are run with a tolerance of eps=0.1.
Denoting eps_lb as the min. epsilon found from the above routines, and eps_ub as the max. epsilon, the following perturbation
radius is chosen: `1/3 eps_lb + 2/3 eps_ub`. 

## How to run the dataset generation

The dataset generation procedure, with a user-specified random seed, can be run via the following command: 
```
# Install the dataset generation code in a new conda environment
source ./installer.sh
# Run the dataset generation
python generate_properties.py --seed <seed>
# Deactivate conda environment
conda deactivate
```
The script generates the .vnnlib properties, the .onnx networks and the .csv containing the benchmark specification.
The dataset generation is executed by default on a GPU (which can be specified by preprending `CUDA_VISIBLE_DEVICES=id`).
In order to instead run the generation on CPU, please add `--cpu` to the command line arguments of `generate_properties.py`. 
 
*NOTE:* On a Nvidia Titan Xp GPU, the complete dataset generation procedure (including the time to download the CIFAR10 test set)
 takes roughly 300s. One can expect slower runtimes on CPU (depending on the number of available cores).


## References.

[1] https://arxiv.org/abs/1711.00851

[2] https://arxiv.org/abs/1710.06081

[3] https://arxiv.org/abs/2011.13824

[4] https://arxiv.org/abs/1811.00866