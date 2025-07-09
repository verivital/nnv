"""Generate benchmarking vnnlib files for trained models."""

import os
import shutil
import argparse
import csv
import random

import torch
import torchvision.transforms as transforms
import torchvision.datasets as datasets

from vnnlib_utils import create_input_bounds, save_vnnlib


models = [
    {
        'name': 'pgd_2_3_16',
        'num_instances': 100,
        'timeout': 100,
    },
    {
        'name': 'ibp_3_3_8',
        'num_instances': 100,
        'timeout': 100,
    }
]

mean = torch.tensor([0.4914, 0.4822, 0.4465])
std = torch.tensor([0.2023, 0.1994, 0.2010])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('seed', type=int)
    args = parser.parse_args()

    random.seed(args.seed)

    for file in ['vnnlib', 'onnx']:
        if os.path.exists(file):
            shutil.rmtree(file)
        os.makedirs(file)

    instances = []
    for model in models:
        name = model['name']
        path = f'models/{name}'
        with open(f'{path}/index.txt') as f:
            indexes = list(map(int, f.readlines()))
        indexes = list(enumerate(indexes))
        random.shuffle(indexes)
        onnx_path = f'onnx/{name}.onnx'
        shutil.copy(f'{path}/model.onnx', onnx_path)
        data = torch.load(f'{path}/data.pt')
        for i in range(model['num_instances']):
            vnnlib_path = f'vnnlib/{name}_{indexes[i][1]}.vnnlib'
            x, y = data[0][indexes[i][0]], data[1][indexes[i][0]]
            input_bounds = create_input_bounds(x, 1./255, mean, std)
            save_vnnlib(input_bounds, y, vnnlib_path, total_output_class=10)
            instances.append((onnx_path, vnnlib_path, model['timeout']))

    with open('instances.csv', 'w') as f:
        csv.writer(f).writerows(instances)
