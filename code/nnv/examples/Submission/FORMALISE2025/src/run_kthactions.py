import os
import sys

import vvn.verify as vvn
from vvn.config import Config


if __name__ == "__main__":
    # Results directory - the prep modules handle subdirectory creation
    # Structure: output_dir/ds_type/ver_algorithm/sample_len/
    output_dir = "/tmp/results/KTHActions"
    os.makedirs(output_dir, exist_ok=True)

    # define the starting configuration
    config = Config(
        class_size=6,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='kthactions',
        sample_len=8,
        ver_algorithm='relax',  # KTH Actions only uses relax algorithm
        timeout=1800,
        output_dir=output_dir
    )

    # 25 samples for KTH Actions
    samples = list(range(1, 26))

    # =====================================
    # KTH Actions Experiments (relax only)
    # =====================================

    print("\n" + "="*50)
    print("Running KTH Actions 8-frame experiments")
    print("="*50 + "\n")
    vvn.run_kthactions(config=config, indices=samples)

    print("\n" + "="*50)
    print("Running KTH Actions 16-frame experiments")
    print("="*50 + "\n")
    config.sample_len = 16
    vvn.run_kthactions(config=config, indices=samples)

    print("\n" + "="*50)
    print("Running KTH Actions 32-frame experiments")
    print("="*50 + "\n")
    config.sample_len = 32
    vvn.run_kthactions(config=config, indices=samples)
