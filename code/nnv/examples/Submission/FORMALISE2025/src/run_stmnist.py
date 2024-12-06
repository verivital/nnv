import os
import sys

import vvn.stmnistprep as vsp
import vvn.verify as vvn
from vvn.config import Config


if __name__ == "__main__":
    # get the results dir
    root = os.path.dirname(os.getcwd())
    output_dir = os.path.join(root, 'FORMALISE2025' ,'results')

    # define the starting configuration 
    config = Config(
        class_size=10,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='stmnist',
        sample_len=16,
        ver_algorithm='relax',
        timeout=1800,
        output_dir=output_dir
    )

    # get the samples you wish to verify
    # samples = vsp.generate_indices(config)
    samples = [311, 64, 719, 646, 595, 385, 268, 220, 1090, 77, 75, 237, 568, 616, 1290, 67, 515, 1086, 579, 1162, 725, 14, 425, 1090, 880, 418, 558, 866, 268, 236, 992, 256, 915, 888, 694, 112, 1181, 1378, 344, 989, 197, 763, 932, 499, 171, 117, 603, 758, 200, 616, 267, 992, 725, 1172, 946, 430, 973, 911, 545, 700, 182, 450, 1373, 646, 431, 1187, 990, 709, 573, 841, 139, 605, 77, 821, 1045, 701, 164, 549, 820, 550, 1282, 1029, 1180, 396, 696, 385, 649, 1382, 693, 1106, 1039, 932, 573, 384, 1298, 1267, 234, 120, 298, 413] 

    if len(sys.argv) < 2:
        samples = samples[:, :, 10]

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

    config = Config(
        class_size=10,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='stmnist',
        sample_len=16,
        ver_algorithm='approx',
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