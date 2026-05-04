import argparse
import os
import sys

import vvn.gtsrbprep as vgp
import vvn.verify as vvn
from vvn.config import Config


def parse_args():
    parser = argparse.ArgumentParser(description='Run GTSRB verification experiments')
    parser.add_argument('--algorithm', type=str, choices=['relax', 'approx', 'both'], default='both',
                        help='Verification algorithm to use (default: both)')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Results directory - the prep modules handle subdirectory creation
    # Structure: output_dir/ds_type/ver_algorithm/sample_len/
    output_dir = "/tmp/results/GTSRB"
    os.makedirs(output_dir, exist_ok=True)

    # get the samples you wish to verify
    samples = list(range(1, 216))

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
            class_size=5,
            epsilon=[1/255, 2/255, 3/255],
            ds_type='gtsrb',
            sample_len=4,
            ver_algorithm=algorithm,
            timeout=1800,
            output_dir=output_dir
        )

        # run experiment #1 : dataset = gtsrb, video length = 4
        vvn.run_gtsrb(config=config, indices=samples)

        # run experiment #2 : dataset = gtsrb, video length = 8
        config.sample_len = 8
        vvn.run_gtsrb(config=config, indices=samples)

        # run experiment #3 : dataset = gtsrb, video length = 16
        config.sample_len = 16
        vvn.run_gtsrb(config=config, indices=samples)
