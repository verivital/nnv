import os
import sys

import vvn.gtsrbprep as vgp
import vvn.verify as vvn
from vvn.config import Config


if __name__ == "__main__":
    # get the results dir
    # root = os.path.dirname(os.getcwd())
    # output_dir = os.path.join(root, 'FORMALISE2025' ,'results')
    # output_dir = '/tmp/GTSRB'
    output_dir = "/tmp/approx_results/GTSRB"

    if not os.path.isdir("/tmp/approx_results/GTSRB"):
        os.mkdir("/tmp/approx_results/GTSRB")
        
    os.mkdir(os.path.join(output_dir, "4"))
    os.mkdir(os.path.join(output_dir, "8"))
    os.mkdir(os.path.join(output_dir, "16"))

    # define the starting configuration 
    config = Config(
        class_size=5,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='gtsrb',
        sample_len=4,
        ver_algorithm='approx',
        timeout=1800,
        output_dir=output_dir
    )

    # get the samples you wish to verify
    samples = list(range(1, 216))

    # if len(sys.argv) < 2:
    #     samples = samples[:, :, 10]

    # =====================================
    # ============ RELAX ==================
    # =====================================

    # run experiment #1 : dataset = gtsrb, video length = 4
    vvn.run_gtsrb(config=config, indices=samples) 

    # run experiment #2 : dataset = gtsrb, video length = 8
    config.sample_len = 8
    vvn.run_gtsrb(config=config, indices=samples)

    # run experiment #3 : dataset = gtsrb, video length = 16
    config.sample_len = 16
    vvn.run_gtsrb(config=config, indices=samples)


    # =====================================
    # ============ APPROX =================
    # =====================================

    # config = Config(
    #     class_size=5,
    #     epsilon=[1/255, 2/255, 3/255],
    #     ds_type='gtsrb',
    #     sample_len=4,
    #     ver_algorithm='approx',
    #     timeout=1800,
    #     output_dir=output_dir
    # )

    # # run experiment #1 : dataset = gtsrb, video length = 4
    # vvn.run_gtsrb(config=config, indices=samples)

    # # run experiment #2 : dataset = gtsrb, video length = 8
    # config.sample_len = 8
    # vvn.run_gtsrb(config=config, indices=samples)

    # # run experiment #3 : dataset = gtsrb, video length = 16
    # config.sample_len = 16
    # vvn.run_gtsrb(config=config, indices=samples)
