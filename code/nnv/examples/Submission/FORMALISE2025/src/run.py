import argparse
import os
import sys

import vvn.prep as vp
import vvn.verify as vvn
from analysis.make_table2 import summarize
from vvn.config import Config


def parse_args():
    parser = argparse.ArgumentParser(description='Run MNIST Video verification experiments')
    parser.add_argument('--algorithm', type=str, choices=['relax', 'approx', 'both'], default='both',
                        help='Verification algorithm to use (default: both)')
    parser.add_argument('--subset', action='store_true',
                        help='Run only a subset of experiments (first row of Table 2)')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    # Results directory - the prep modules handle subdirectory creation
    # Structure: output_dir/ds_type/ver_algorithm/sample_len/
    output_dir = "/tmp/results/MNIST"
    os.makedirs(output_dir, exist_ok=True)

    # get the samples you wish to verify
    zoom_in_samples = list(range(1, 101))
    zoom_out_samples = list(range(1, 101))

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
            ds_type='zoom_in',
            sample_len=4,
            ver_algorithm=algorithm,
            timeout=1800,
            output_dir=output_dir
        )

        # run experiment #1 : dataset = zoom in, video length = 4
        vvn.run(config=config, indices=zoom_in_samples)

        if args.subset:
            # summarize the results
            tv, rt = summarize(vp.build_output_filepath(config=config, parent_only=True))
            rt = [round(rt[i], 2) for i in range(len(rt))]

            label = f"Zoom In, {algorithm}     1/255       2/255       3/255       "
            tv_row = f"{'PSRV':>{14}}" + f"{' ':<{5}}" + f"{tv[0]:<{12}}" + f"{tv[1]:<{12}}" + f"{tv[2]:<{12}}"
            rt_row = f"{'Avg. Time (s)':>{14}}" + f"{' ':<{5}}" + f"{rt[0]:<{12}}" + f"{rt[1]:<{12}}" + f"{rt[2]:<{12}}"

            length_toplines = max(len(tv_row), len(rt_row), len(label))

            # print results
            print('\n')
            print("="*length_toplines)
            print(label)
            print("-"*length_toplines)
            print(tv_row)
            print(rt_row)
            print("="*length_toplines)
            print('\n')

            # move on to next algorithm if running both
            continue

        # run experiment #2 : dataset = zoom out, video length = 4
        config.ds_type = 'zoom_out'
        vvn.run(config=config, indices=zoom_out_samples)

        # run experiment #3 : dataset = zoom in , video length = 8
        config.ds_type = 'zoom_in'
        config.sample_len = 8
        vvn.run(config=config, indices=zoom_in_samples)

        # run experiment #4 : dataset = zoom out, video length = 8
        config.ds_type = 'zoom_out'
        vvn.run(config=config, indices=zoom_out_samples)

        # run experiment #5 : dataset = zoom in, video length = 16
        config.ds_type = 'zoom_in'
        config.sample_len = 16
        vvn.run(config=config, indices=zoom_in_samples)

        # run experiment #6 : dataset = zoom out, video length = 16
        config.ds_type = 'zoom_out'
        vvn.run(config=config, indices=zoom_out_samples)
