import os
import sys

import vvn.prep as vp
import vvn.verify as vvn
from vvn.config import Config


if __name__ == "__main__":
    # get the results dir
    root = os.path.dirname(os.getcwd()) # Submission folder
    output_dir = os.path.join(root, 'FORMALISE2025' ,'results')

    # define the starting configuration 
    config = Config(
        class_size=10,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='zoom_in',
        sample_len=4,
        ver_algorithm='relax',
        timeout=1800,
        output_dir=output_dir
    )

    # get the samples you wish to verify
    # zoom_in_samples, zoom_out_samples = vp.generate_indices(config)
    zoom_in_samples = [751, 137, 34, 878, 328, 289, 262, 168, 871, 126, 796, 877, 645, 108, 693, 501, 44, 41, 117, 257, 272, 597, 707, 37, 661, 237, 845, 763, 830, 645, 498, 259, 532, 692, 331, 972, 8, 904, 966, 193, 826, 501, 407, 331, 187, 254, 908, 402, 126, 116, 452, 121, 429, 411, 709, 317, 969, 57, 863, 543, 636, 150, 450, 99, 652, 350, 999, 736, 724, 432, 680, 230, 833, 90, 61, 780, 267, 922, 346, 100, 272, 125, 452, 331, 537, 744, 435, 198, 442, 423, 248, 790, 320, 830, 806, 761, 92, 714, 744, 207] 
    zoom_out_samples = [624, 882, 278, 180, 540, 439, 306, 757, 821, 654, 248, 817, 368, 949, 963, 59, 260, 34, 357, 465, 304, 69, 238, 666, 867, 356, 239, 776, 585, 460, 760, 536, 158, 301, 154, 280, 908, 659, 632, 297, 910, 687, 499, 686, 463, 418, 248, 152, 596, 578, 96, 922, 50, 117, 169, 738, 176, 989, 809, 491, 702, 67, 445, 441, 547, 616, 285, 649, 12, 809, 872, 126, 812, 630, 916, 303, 952, 758, 390, 120, 332, 507, 174, 529, 4, 873, 868, 297, 586, 933, 196, 594, 112, 736, 337, 755, 719, 223, 169, 433] 

    if len(sys.argv) < 2:
        zoom_in_samples = zoom_in_samples[:, :, 10]
        zoom_out_samples = zoom_out_samples[:, :, 10]

    # =====================================
    # ============ RELAX ==================
    # =====================================

    # run experiment #1 : dataset = zoom in, video length = 4
    vvn.run(config=config, indices=zoom_in_samples)

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


    # =====================================
    # ============ APPROX =================
    # =====================================

    config = Config(
        class_size=10,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='zoom_in',
        sample_len=4,
        ver_algorithm='approx',
        timeout=1800,
        output_dir=output_dir
    )

    # run experiment #1 : dataset = zoom in, video length = 4
    vvn.run(config=config, indices=zoom_in_samples)

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

