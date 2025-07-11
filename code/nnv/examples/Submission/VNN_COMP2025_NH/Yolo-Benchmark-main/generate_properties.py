import os
import argparse
import csv
import random
import sys

import numpy as np
import torch
import torch.nn.functional as F
import torchvision.datasets as dset
import torchvision.transforms as trans
from torch.utils.data import DataLoader
from torch.utils.data import sampler
from gen_vnnlib import gen_vnnlib

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('seed', type=int, default=0, help='random seed.')
    args = parser.parse_args()
    args.selected_instances = 72
    args.timeout = 300

    input_path = 'inputs'
    output_path = 'raw_outputs'
    model_name = 'TinyYOLO.onnx'

    assert os.path.exists(os.path.join('onnx', model_name))

    if not os.path.exists('vnnlib'):
        os.makedirs('vnnlib')
    for vnnlib_file in os.scandir('vnnlib'):
        os.remove(vnnlib_file.path)
    if os.path.exists('instances.csv'):
        os.remove('instances.csv')

    sample_list = os.listdir(input_path)

    np.random.seed(args.seed)
    selected = np.random.permutation(len(sample_list))[:args.selected_instances]
    with open('instances.csv', 'w', newline='') as f:
        csv_writer = csv.writer(f)
        for selected_ind in selected:
            im_name = sample_list[selected_ind].split('.')[0]
            im = torch.load(os.path.join(input_path, sample_list[selected_ind]))
            # print(os.path.join(output_path, sample_list[selected_ind]))
            raw_output = torch.load(os.path.join(output_path, sample_list[selected_ind]))
            gen_vnnlib(im, raw_output, im_name)
            csv_writer.writerow(['onnx/' + model_name, 'vnnlib/TinyYOLO_prop_{}_eps_1_255.vnnlib'.format(im_name), args.timeout])
            