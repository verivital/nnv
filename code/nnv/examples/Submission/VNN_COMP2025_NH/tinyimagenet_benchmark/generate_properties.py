import os
import argparse
import csv
import random
import sys

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data import sampler

tinyimagenet_mean = (0.4802, 0.4481, 0.3975)
tinyimagenet_std = (0.2302, 0.2265, 0.2262)


# noinspection PyShadowingNames
def create_input_bounds(img: torch.Tensor, eps: float,
                        mean: tuple, std: tuple) -> torch.Tensor:
    """
    Creates input bounds for the given image and epsilon.

    The lower bounds are calculated as img-eps clipped to [0, 1] and the upper bounds
    as img+eps clipped to [0, 1].

    Args:
        img:
            The image.
        eps:
           The maximum accepted epsilon perturbation of each pixel.
        mean:
            The channel-wise means.
        std:
            The channel-wise standard deviation.
    Returns:
        A  img.shape x 2 tensor with the lower bounds in [..., 0] and upper bounds
        in [..., 1].
    """

    mean = torch.tensor(mean, device=img.device).view(-1, 1, 1)
    std = torch.tensor(std, device=img.device).view(-1, 1, 1)

    bounds = torch.zeros((*img.shape, 2), dtype=torch.float32)
    bounds[..., 0] = (np.clip((img - eps), 0, 1) - mean) / std
    bounds[..., 1] = (np.clip((img + eps), 0, 1) - mean) / std
    # print(bounds[..., 0].abs().sum(), bounds[..., 1].abs().sum())

    return bounds.view(-1, 2)


# noinspection PyShadowingNames
def save_vnnlib(input_bounds: torch.Tensor, label: int, spec_path: str, total_output_class: int):
    """
    Saves the classification property derived as vnn_lib format.

    Args:
        input_bounds:
            A Nx2 tensor with lower bounds in the first column and upper bounds
            in the second.
        label:
            The correct classification class.
        spec_path:
            The path used for saving the vnn-lib file.
        total_output_class:
            The total number of classification classes.
    """

    with open(spec_path, "w") as f:

        f.write(f"; Property with label: {label}.\n")

        # Declare input variables.
        f.write("\n")
        for i in range(input_bounds.shape[0]):
            f.write(f"(declare-const X_{i} Real)\n")
        f.write("\n")

        # Declare output variables.
        f.write("\n")
        for i in range(total_output_class):
            f.write(f"(declare-const Y_{i} Real)\n")
        f.write("\n")

        # Define input constraints.
        f.write(f"; Input constraints:\n")
        for i in range(input_bounds.shape[0]):
            f.write(f"(assert (<= X_{i} {input_bounds[i, 1]}))\n")
            f.write(f"(assert (>= X_{i} {input_bounds[i, 0]}))\n")
            f.write("\n")
        f.write("\n")

        # Define output constraints.
        f.write(f"; Output constraints:\n")
        # orignal separate version:
        # for i in range(total_output_class):
        #     if i != label:
        #         f.write(f"(assert (>= Y_{label} Y_{i}))\n")
        # f.write("\n")

        # disjunction version:
        f.write("(assert (or\n")
        for i in range(total_output_class):
            if i != label:
                f.write(f"    (and (>= Y_{i} Y_{label}))\n")
        f.write("))\n")



def create_vnnlib(args, dataset):

    assert os.path.exists(f"./onnx/{args.dataset}_resnet_{args.model}.onnx")
    instance_list = []
    epsilons = [eval(eps) for eps in args.epsilons.split(" ")]

    init_dir = f"./{args.mode}/".replace("generate", "generated").replace("_csv", "")
    if not os.path.isdir(init_dir):
        os.mkdir(init_dir)

    result_dir = init_dir

    print(result_dir)

    if not os.path.isdir(result_dir):
        os.mkdir(result_dir)

    mu = torch.tensor(tinyimagenet_mean).view(3, 1, 1)
    std = torch.tensor(tinyimagenet_std).view(3, 1, 1)

    normalize = lambda X: (X - mu) / std

    np.random.seed(args.seed)
    random.shuffle(dataset)
    dataset = dataset[:args.selected_instances]

    for eps in epsilons:
        cnt = 0
        for image, label, sidx, idx in dataset:
            image = image.unsqueeze(0)
            y = torch.tensor([label], dtype=torch.int64)

            print("scanned images: {}, selected: {}, label {}".format(sidx, cnt, label))

            input_bounds = create_input_bounds(image, eps, mean = mu, std = std)
            vnnlib_path = f"{args.dataset}_resnet_{args.model}_prop_idx_{idx}_sidx_{sidx}_eps_{eps:.4f}.vnnlib"
            spec_path = os.path.join(result_dir, vnnlib_path)
            save_vnnlib(input_bounds, label, spec_path, total_output_class=200)

            instance_list.append([
                f"{args.dataset}_resnet_{args.model}.onnx",
                vnnlib_path, f"{args.timeout}"])
            cnt += 1

    assert os.path.exists(f"./generated_vnnlib/")

    with open(f'./instances.csv', 'a') as f:
        write = csv.writer(f)
        write.writerows(instance_list)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('seed', type=int, default=0, help='random seed.') # seed for points selection
    args = parser.parse_args()
    args.epsilons = "1/255"
    args.mode = "generate_vnnlib_csv"
    vnnlib_path = "generated_vnnlib"
    if not os.path.exists(vnnlib_path):
        os.makedirs(vnnlib_path)
    # Remove old files in the vnnlib folder.
    for vnnlib_file in os.scandir(vnnlib_path):
        os.remove(vnnlib_file.path)
    if os.path.exists("instances.csv"):
        os.remove("instances.csv")
    dataset = "TinyImageNet"
    model = "medium"
    args.dataset = "TinyImageNet"
    args.model = model
    args.selected_instances = 200
    args.timeout = 100
    dataset = torch.load(f'data/{args.dataset}_resnet_{args.model}.pt')
    create_vnnlib(args, dataset)
