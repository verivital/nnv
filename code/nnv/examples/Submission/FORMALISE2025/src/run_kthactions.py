import os
import sys

import vvn.verify as vvn
from vvn.config import Config

if __name__=="__main__":
    output_dir = "/tmp"

    # define the starting configurations
    config = Config(
        class_size=6,
        epsilon=[1/255,2/255,3/255],
        ds_type='kthactions',
        sample_len=8,
        ver_algorithm='relax',
        timeout=1800,
        output_dir=output_dir
    )

    samples = list(range(1, 26))
    vvn.run_kthactions(config=config, indices=samples)

    config.sample_len = 16
    vvn.run_kthactions(config=config, indices=samples)

    config.sample_len = 32
    vvn.run_kthactions(config=config, indices=samples)
