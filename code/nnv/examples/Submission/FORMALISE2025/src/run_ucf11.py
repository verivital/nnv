import os
import sys

import vvn.verify as vvn
from vvn.config import Config


if __name__ == "__main__":
    # Results directory - the prep modules handle subdirectory creation
    # Structure: output_dir/ds_type/ver_algorithm/sample_len/
    output_dir = "/tmp/results/UCF11"
    os.makedirs(output_dir, exist_ok=True)

    # define the starting configuration
    config = Config(
        class_size=11,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='ucf11',
        sample_len=8,
        ver_algorithm='relax',  # UCF11 only uses relax algorithm
        timeout=1800,
        output_dir=output_dir
    )

    # 25 samples for UCF11
    samples = list(range(1, 26))

    # =====================================
    # Table 5: UCF11 Experiments (relax only)
    # =====================================

    #
    # 8 frames - all output channel sizes
    #
    print("\n" + "="*50)
    print("Running UCF11 8-frame experiments")
    print("="*50 + "\n")
    vvn.run_ucf11(config=config, indices=samples, out_channels=8)
    vvn.run_ucf11(config=config, indices=samples, out_channels=16)
    vvn.run_ucf11(config=config, indices=samples, out_channels=32)
    vvn.run_ucf11(config=config, indices=samples, out_channels=64)

    #
    # 16 frames - output channels 8, 16, 32
    #
    print("\n" + "="*50)
    print("Running UCF11 16-frame experiments")
    print("="*50 + "\n")
    config.sample_len = 16
    vvn.run_ucf11(config=config, indices=samples, out_channels=8)
    vvn.run_ucf11(config=config, indices=samples, out_channels=16)
    vvn.run_ucf11(config=config, indices=samples, out_channels=32)

    #
    # 32 frames - output channels 8, 16
    #
    print("\n" + "="*50)
    print("Running UCF11 32-frame experiments")
    print("="*50 + "\n")
    config.sample_len = 32
    vvn.run_ucf11(config=config, indices=samples, out_channels=8)
    vvn.run_ucf11(config=config, indices=samples, out_channels=16)
