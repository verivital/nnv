# CGAN benchmark for VNNCOMP 2023

Run
----------------------

```bash
python generate_instances.py <seed>
```

Installation and requirement
----------------------

```bash
pip install -r requirements.txt
```


Description
----------------------
This is a benchmark of conditional generative adversarial networks (cGANs) for VNNCOMP 2023.

The objective of this cGAN is to generate camera images that contain a vehicle obstacle located at a specific distance in front of the ego vehicle, where the distance is controlled by the input distance condition.

Here we attach some generated images as well as the architecture of the cGAN (including both generator and discriminator).

<p align="center">
    <img alt="output_images" src="https://github.com/stanleybak/vnncomp2023/assets/29678149/fd6a36fd-42dd-4361-b28f-c485d494c87e">
</p>
<p align="center">    
    <img alt="cgan_architecture" src="https://github.com/stanleybak/vnncomp2023/assets/29678149/51a653bd-dbdd-4944-ab14-23be3370f565">
</p>

The generator takes two inputs: 1) a distance condition (1-d scalar) and 2) a noise vector that controls the environment (4-d vector). The output of the generator is a generated image.

The discriminator takes the generated image as input and outputs two values: 1) a real/fake score (1-d scalar) and 2) a predicted distance (1-d scalar).

For verification, we combine these two components together and set proper verification specifications for input distance, input noise, and predicted distance.

We also offer several different models with varying architectures (CNN and Transformer) and image sizes (32x32, 64x64) to provide a range of difficulty levels.
