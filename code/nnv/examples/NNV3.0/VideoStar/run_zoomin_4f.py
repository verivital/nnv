#!/usr/bin/env python3
"""
VideoStar ZoomIn-4f Verification Script

Optional Python entry point that drives the same MATLAB verification as
``run_zoomin_4f.m`` via the MATLAB engine for Python. The .m file is the
canonical reference implementation; this wrapper is provided for users
who want a Python-native CLI.

Usage:
    python run_zoomin_4f.py [--algorithm {relax,approx,both}] [--num-samples N]

Requirements:
    - MATLAB R2024b with the Python engine installed (`matlab.engine`)
    - NNV toolbox available at <repo>/code/nnv (auto-resolved below)
    - All assets bundled in this folder (models/, data/ZoomIn/, src/vvn/, npy-matlab/)
"""

import argparse
import csv
import datetime
import io
import os
import sys
from typing import Tuple

import matlab.engine


# Paths resolved relative to this file so the script is self-contained.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
NNV_PATH = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', '..'))
VVN_SRC_PATH = os.path.join(SCRIPT_DIR, 'src', 'vvn')
NPY_MATLAB_PATH = os.path.join(SCRIPT_DIR, 'npy-matlab')
RESULTS_DIR = os.path.join(
    SCRIPT_DIR, 'results', datetime.datetime.now().strftime('%y%m%d-%H%M%S')
)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Run ZoomIn-4f verification experiments for VideoStar'
    )
    parser.add_argument(
        '--algorithm',
        type=str,
        choices=['relax', 'approx', 'both'],
        default='relax',
        help='Verification algorithm to use (default: relax)'
    )
    parser.add_argument(
        '--num-samples',
        type=int,
        default=10,
        help='Number of samples to verify (default: 10)'
    )
    parser.add_argument(
        '--timeout',
        type=int,
        default=1800,
        help='Timeout per sample in seconds (default: 1800)'
    )
    return parser.parse_args()


def prepare_engine(nnv_path: str, vvn_src_path: str, npy_matlab_path: str):
    """Start MATLAB engine and add required paths."""
    eng = matlab.engine.start_matlab()
    print('Started MATLAB engine!')

    # Add paths (NNV toolbox + vendored vvn + npy-matlab).
    eng.addpath(eng.genpath(nnv_path))
    eng.addpath(vvn_src_path)
    eng.addpath(npy_matlab_path)

    # verifyvideo.m uses relative paths (`data/ZoomIn/...`, `models/...`),
    # so cd into VideoStar/ where those exist.
    eng.cd(SCRIPT_DIR)

    return eng


def verify(
    ds_type: str,
    sample_len: int,
    ver_algorithm: str,
    eng,
    index: int,
    eps_index: int,
    timeout: int
) -> Tuple[int, float, str]:
    """Run verification on a single sample."""
    if not eng:
        raise Exception(
            'MATLAB Engine was not correctly started. '
            'Please make sure to run `prepare_engine`.'
        )

    # Call MATLAB script to run verification
    future = eng.verifyvideo(
        ds_type,
        sample_len,
        ver_algorithm,
        index,
        eps_index,
        nargout=3,
        background=True,
        stdout=io.StringIO()
    )

    try:
        [res, t, met] = future.result(timeout=float(timeout))
    except matlab.engine.TimeoutError:
        print('  Timeout')
        res = 3
        t = 'timeout'
        met = 'timeout'

    future.cancel()

    return res, t, met


def write_results(output_file: str, sample_num: int, res: int, t, met: str):
    """Append results to CSV file."""
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write header if file doesn't exist
    if not os.path.exists(output_file):
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Sample Number', 'Result', 'Time', 'Method'])

    # Append results
    with open(output_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([sample_num, res, t, met])


def run_zoomin_4f(algorithm: str, num_samples: int, timeout: int):
    """Run ZoomIn-4f verification experiments."""
    # Configuration
    ds_type = 'zoom_in'
    sample_len = 4
    epsilon_indices = [1, 2, 3]  # 1/255, 2/255, 3/255

    # Sample indices to verify (1-indexed for MATLAB)
    sample_indices = list(range(1, num_samples + 1))

    # Output directory: results/<timestamp>/
    output_dir = RESULTS_DIR
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 50)
    print(f"Running ZoomIn-4f verification")
    print("=" * 50)
    print(f"  Algorithm: {algorithm}")
    print(f"  Dataset: {ds_type}")
    print(f"  Sample length: {sample_len} frames")
    print(f"  Number of samples: {num_samples}")
    print(f"  Timeout: {timeout} seconds")
    print(f"  Output directory: {output_dir}")
    print()

    # Start MATLAB engine
    eng = prepare_engine(NNV_PATH, VVN_SRC_PATH, NPY_MATLAB_PATH)

    # Run verification
    for sample_idx, index in enumerate(sample_indices):
        print(f'Sample {index} ({sample_idx + 1}/{num_samples})')

        for eps_index in epsilon_indices:
            output_file = os.path.join(output_dir, f'eps={eps_index}_255.csv')

            # Verify the sample
            res, t, met = verify(
                ds_type, sample_len, algorithm, eng, index, eps_index, timeout
            )

            # Write results
            write_results(output_file, index, res, t, met)

            # Print result
            if isinstance(t, (int, float)):
                print(f'  eps={eps_index}/255: Result={res}, Time={t:.2f}s')
            else:
                print(f'  eps={eps_index}/255: Result={res}, Time={t}')

    # Close MATLAB engine
    eng.quit()

    print()
    print("=" * 50)
    print("Verification complete.")
    print(f"Results saved to: {output_dir}")
    print("=" * 50)


def main():
    args = parse_args()

    # Determine which algorithms to run
    algorithms = []
    if args.algorithm in ['relax', 'both']:
        algorithms.append('relax')
    if args.algorithm in ['approx', 'both']:
        algorithms.append('approx')

    for algorithm in algorithms:
        run_zoomin_4f(algorithm, args.num_samples, args.timeout)


if __name__ == "__main__":
    main()
