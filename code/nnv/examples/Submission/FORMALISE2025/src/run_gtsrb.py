import os
import sys

import vvn.gtsrbprep as vgp
import vvn.verify as vvn
from vvn.config import Config


if __name__ == "__main__":
    # get the results dir
    root = os.path.dirname(os.getcwd())
    output_dir = os.path.join(root, 'FORMALISE2025' ,'results')

    # define the starting configuration 
    config = Config(
        class_size=5,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='gtsrb',
        sample_len=4,
        ver_algorithm='relax',
        timeout=1800,
        output_dir=output_dir
    )

    # get the samples you wish to verify
    # samples = vgp.generate_indices(config)
    samples = [12597, 2196, 507, 5423, 4840, 4410, 2743, 2020, 10754, 1707, 11622, 8312, 635, 596, 1851, 4320, 4595, 9948, 11849, 535, 11062, 3942, 10746, 8256, 4356, 8857, 11599, 5483, 127, 3154, 8326, 6696, 5474, 3075, 4256, 6625, 2013, 1832, 7475, 1908, 7068, 6765, 11879, 5214, 858, 9051, 10579, 2459, 7450, 1550, 10884, 5771, 12391, 12180, 7128, 11366, 3809, 1367, 910, 4497, 5706, 1570, 4599, 1989, 7478, 5480, 8937, 12524, 7190, 3209, 7303, 6996, 4143, 5263, 1406, 11985, 12510, 3378, 10534, 4838, 3228, 9106, 7467, 5324, 12607, 10979, 4343, 6387, 1109, 4525, 641, 6208, 7898, 5278, 1302, 4171, 11173, 6190, 4204, 9824, 7788, 9041, 2814, 5221, 2744, 4873, 11065, 10632, 5180, 11510, 8452, 11492, 7863, 7133, 4337, 2721, 10042, 9713, 1796, 935, 2161, 3020, 12365, 3161, 8318, 11739, 1249, 7568, 7504, 11728, 9222, 10436, 4971, 10910, 228, 2253, 10594, 5260, 12629, 6695, 2199, 5778, 8578, 3125, 8944, 61, 5195, 9852, 3527, 10001, 2101, 12318, 5875, 12588, 10001, 11982, 3943, 3020, 7371, 3190, 10639, 10457, 10, 11787, 6384, 9622, 398, 2204, 7152, 6052, 4732, 1147, 4760, 11178, 1549, 1685, 9577, 1358, 10501, 2476, 2531, 9366, 10839, 3267, 5224, 10405, 11933, 8336, 4190, 10637, 3985, 6135, 7857, 7363, 8642, 10196, 8901, 2385, 4901, 4439, 1258, 6656, 427, 11580, 10921, 4549, 11583, 4353, 141, 1398, 12441, 1162, 4523, 1326, 629] 

    if len(sys.argv) < 2:
        samples = samples[:, :, 10]

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

    config = Config(
        class_size=5,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='gtsrb',
        sample_len=4,
        ver_algorithm='approx',
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
