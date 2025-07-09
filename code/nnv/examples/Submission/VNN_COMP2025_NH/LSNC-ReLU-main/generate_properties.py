"""Generate benchmarking vnnlib files for trained models."""

import os
import shutil
import argparse
import csv
import random
import torch


tolerance = 1e-6

models = [
    {
        'name': 'quadrotor2d_state',
        'num_instances': 80,
        'timeout': 25,
        'box_limit': torch.tensor([0.7, 0.7, 1, 3, 3, 2]),
        'eps': 0.05,
    }
]


def generate_instance(model, vnnlib_path):
    if random.random() < 0.2:
        c1 = random.uniform(0.475, 0.95)
    else:
        c1 = random.uniform(0.05, 0.425)
    c2 = c1 + model['eps']
    generate_vnnlib(model, -model['box_limit'], model['box_limit'], c1, c2, vnnlib_path)
    
def generate_complete_instance(c1, model, vnnlib_path):
    c2 = c1 + model['eps']
    generate_vnnlib(model, -model['box_limit'], model['box_limit'], c1, c2, vnnlib_path)


def generate_vnnlib(model, lower, upper, c1, c2, vnnlib_path):
    state_dim = len(lower)
    with open(vnnlib_path, 'w') as out:
        for i in range(state_dim):
            out.write(f"(declare-const X_{i} Real)\n")
        out.write("(declare-const Y_0 Real)\n")
        out.write("(declare-const Y_1 Real)\n")
        for i in range(2, 2 + state_dim):
            out.write(f"(declare-const Y_{i} Real)\n")
        out.write("\n")
        for i, (l, u) in enumerate(zip(lower, upper)):
            out.write(f"(assert (<= X_{i} {u}))\n")
            out.write(f"(assert (>= X_{i} {l}))\n\n")
        out.write(f"(assert (or\n")
        out.write(f"  (and (>= Y_0 {tolerance}))\n")
        for i in range(state_dim):
            out.write(f"  (and (<= Y_{i+2} {-model['box_limit'][i] - tolerance}))\n")
            out.write(f"  (and (>= Y_{i+2} {model['box_limit'][i] + tolerance}))\n")
        out.write("))\n")
        out.write(f'(assert (<= Y_1 {c2}))\n')
        out.write(f'(assert (>= Y_1 {c1}))\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('seed', type=int)
    args = parser.parse_args()

    random.seed(args.seed)

    if os.path.exists('vnnlib'):
        shutil.rmtree('vnnlib')
    os.makedirs('vnnlib')

    instances = []
    for model in models:
        for i in range(model['num_instances']):
            vnnlib_path = f'vnnlib/{model["name"]}_{i}.vnnlib'
            generate_instance(model, vnnlib_path)
            instances.append(vnnlib_path)

    model_path = 'onnx/relu_quadrotor2d_state.onnx'

    with open('instances.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows([[model_path, path, model["timeout"]] for path in instances])
