## TinyImageNet Neural Network Verification Benchmark for VNN-COMP 2024

This is the re-proposed TinyImageNet benchmark for VNN-COMP 2024. The original benchmark [cifar100_tinyimagenet](https://github.com/Lucas110550/CIFAR100_TinyImageNet_ResNet) was from VNN-COMP 2022.

To reduce complexity to improve its accessibility to new tools and new VNN-COMP participants,
this TinyImageNet model is separated as a different benchmark, since not all tools can handle large model/data like this.
The model contain FC, conv, and ReLU layers only.

For more information about the models, you can check out the original repo in 2022 (here, we only kept the `TinyImageNet_resnet_medium` model only; other models were removed to reduce complexity).

### Setup 

PyTorch is required, which can be installed by `pip install torch`.

To generate specficiations:
```bash
python generate_properties.py SEED
```
