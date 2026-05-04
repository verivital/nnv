import argparse
import os
import sys

import vvn.stmnistprep as vsp
import vvn.verify as vvn
from vvn.config import Config


def parse_args():
    parser = argparse.ArgumentParser(description='Run ST-MNIST verification experiments')
    parser.add_argument('--algorithm', type=str, choices=['relax', 'approx', 'both'], default='both',
                        help='Verification algorithm to use (default: both)')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Results directory - the prep modules handle subdirectory creation
    # Structure: output_dir/ds_type/ver_algorithm/sample_len/
    output_dir = '/tmp/results/STMNIST'
    os.makedirs(output_dir, exist_ok=True)

    # get the samples you wish to verify
    samples = list(range(1, 101))

    # Determine which algorithms to run
    algorithms = []
    if args.algorithm in ['relax', 'both']:
        algorithms.append('relax')
    if args.algorithm in ['approx', 'both']:
        algorithms.append('approx')

    for algorithm in algorithms:
        print(f"\n{'='*50}")
        print(f"Running {algorithm.upper()} algorithm experiments")
        print(f"{'='*50}\n")

        # define the starting configuration
        config = Config(
            class_size=10,
            epsilon=[1/255, 2/255, 3/255],
            ds_type='stmnist',
            sample_len=16,
            ver_algorithm=algorithm,
            timeout=1800,
            output_dir=output_dir
        )

        # run experiment #1 : dataset = stmnist, video length = 16
        vvn.run_stmnist(config=config, indices=samples)

        # run experiment #2 : dataset = stmnist, video length = 32
        config.sample_len = 32
        vvn.run_stmnist(config=config, indices=samples)

        # run experiment #3 : dataset = stmnist, video length = 64
        config.sample_len = 64
        vvn.run_stmnist(config=config, indices=samples)
