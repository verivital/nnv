import os
import sys

import vvn.stmnistprep as vsp
import vvn.verify as vvn
from vvn.config import Config


if __name__ == "__main__":
    # get the results dir
    # root = os.path.dirname(os.getcwd())
    # output_dir = os.path.join(root, 'FORMALISE2025' ,'results')
    # output_dir = '/tmp/STMNIST'
    output_dir = '/tmp/approx_results/STMNIST'

    if not os.path.isdir("/tmp/approx_results/STMNIST"):
        os.mkdir("/tmp/approx_results/STMNIST")

    os.mkdir(os.path.join(output_dir, "16"))
    os.mkdir(os.path.join(output_dir, "32"))
    os.mkdir(os.path.join(output_dir, "64"))

    # define the starting configuration 
    config = Config(
        class_size=10,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='stmnist',
        sample_len=16,
        ver_algorithm='approx',
        timeout=1800,
        output_dir=output_dir
    )

    # get the samples you wish to verify
    samples = list(range(1, 101))

    # if len(sys.argv) < 2:
    #     samples = samples[:, :, 10]

    # =====================================
    # ============ RELAX ==================
    # =====================================

    # run experiment #1 : dataset = stmnist, video length = 16
    vvn.run_stmnist(config=config, indices=samples)

    # run experiment #2 : dataset = stmnist, video length = 32
    config.sample_len = 32
    vvn.run_stmnist(config=config, indices=samples)

    # run experiment #3 : dataset = stmnist, video length = 64
    config.sample_len = 64
    vvn.run_stmnist(config=config, indices=samples)

    # =====================================
    # ============ APPROX =================
    # =====================================

    # config = Config(
    #     class_size=10,
    #     epsilon=[1/255, 2/255, 3/255],
    #     ds_type='stmnist',
    #     sample_len=16,
    #     ver_algorithm='approx',
    #     timeout=1800,
    #     output_dir=output_dir
    # )

    # # run experiment #1 : dataset = stmnist, video length = 16
    # vvn.run_stmnist(config=config, indices=samples)

    # # run experiment #2 : dataset = stmnist, video length = 32
    # config.sample_len = 32
    # vvn.run_stmnist(config=config, indices=samples)

    # # run experiment #3 : dataset = stmnist, video length = 64
    # config.sample_len = 64
    # vvn.run_stmnist(config=config, indices=samples)