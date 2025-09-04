import os
import sys

import vvn.verify as vvn
from vvn.config import Config

if __name__=="__main__":
    # get the results dir
    root = os.path.dirname(os.getcwd())
    output_dir = os.path.join(root, 'FORMALISE2025', 'results')

    # define the starting configurations
    config = Config(
        class_size=11,
        epsilon=[1/255,2/255,3/255],
        ds_type='ucf11',
        sample_len=16,
        ver_algorithm='relax',
        timeout=1800,
        output_dir=output_dir
    )

    # samples = list(range(1, 101))
    samples = [1]
    vvn.run_ucf11(config=config, indices=samples)

    # config.sample_len = 16
    # vvn.run_ucf11(config=config, indices=samples)

    # config.sample_len = 32
    # vvn.run_ucf11(config=config, indices=samples)
