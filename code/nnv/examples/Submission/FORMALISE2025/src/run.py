import os
import sys

import vvn.prep as vp
import vvn.verify as vvn
from analysis.make_table2 import summarize
from vvn.config import Config


if __name__ == "__main__":
    # get the results dir
    # root = os.path.dirname(os.getcwd()) # Submission folder
    # output_dir = os.path.join(root, 'FORMALISE2025' ,'results')
    # output_dir = "/tmp/MNIST"
    output_dir = "/tmp/approx_results/MNIST"

    if not os.path.isdir("/tmp/approx_results/MNIST"):
        os.mkdir("/tmp/approx_results/MNIST")

    # make sure to make these are not already created
    
    # split between the datasets
    os.mkdir(os.path.join(output_dir, "ZoomIn"))
    os.mkdir(os.path.join(output_dir, "ZoomOut"))
    
    # split between number of frames
    os.mkdir(os.path.join(output_dir, "ZoomIn", "4"))
    os.mkdir(os.path.join(output_dir, "ZoomIn", "8"))
    os.mkdir(os.path.join(output_dir, "ZoomIn", "16"))

    os.mkdir(os.path.join(output_dir, "ZoomOut", "4"))
    os.mkdir(os.path.join(output_dir, "ZoomOut", "8"))
    os.mkdir(os.path.join(output_dir, "ZoomOut", "16"))

    # define the starting configuration 
    config = Config(
        class_size=10,
        epsilon=[1/255, 2/255, 3/255],
        ds_type='zoom_in',
        sample_len=4,
        ver_algorithm='approx',
        timeout=1800,
        output_dir=output_dir
    )

    # get the samples you wish to verify  
    zoom_in_samples = list(range(1, 101))
    zoom_out_samples = list(range(1, 101))

    # =====================================
    # ============ RELAX ==================
    # =====================================

    # run experiment #1 : dataset = zoom in, video length = 4
    vvn.run(config=config, indices=zoom_in_samples)

    if "--subset" in sys.argv:
        # summarize the results
        tv, rt = summarize(vp.build_output_filepath(config=config, parent_only=True))
        rt = [round(rt[i], 2) for i in range(len(rt))]

        label = "Zoom In, relax     1/255       2/255       3/255       "
        tv_row = f"{'PSRV':>{14}}" + f"{' ':<{5}}" + f"{tv[0]:<{12}}" + f"{tv[1]:<{12}}" + f"{tv[2]:<{12}}"
        rt_row = f"{'Avg. Time (s)':>{14}}" + f"{' ':<{5}}" + f"{rt[0]:<{12}}" + f"{rt[1]:<{12}}" + f"{rt[2]:<{12}}"

        length_toplines = max(len(tv_row), len(rt_row), len(label))

        # print results
        print('\n')
        print("="*length_toplines)
        print(label)
        print("-"*length_toplines)
        print(tv_row)
        print(rt_row)
        print("="*length_toplines)
        print('\n')

        # ensure the program stops there
        sys.exit()

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

    # config = Config(
    #     class_size=10,
    #     epsilon=[1/255, 2/255, 3/255],
    #     ds_type='zoom_in',
    #     sample_len=4,
    #     ver_algorithm='approx',
    #     timeout=1800,
    #     output_dir=output_dir
    # )

    # # run experiment #1 : dataset = zoom in, video length = 4
    # vvn.run(config=config, indices=zoom_in_samples)

    # # run experiment #2 : dataset = zoom out, video length = 4
    # config.ds_type = 'zoom_out'
    # vvn.run(config=config, indices=zoom_out_samples)

    # # run experiment #3 : dataset = zoom in , video length = 8
    # config.ds_type = 'zoom_in'
    # config.sample_len = 8
    # vvn.run(config=config, indices=zoom_in_samples)

    # # run experiment #4 : dataset = zoom out, video length = 8
    # config.ds_type = 'zoom_out'
    # vvn.run(config=config, indices=zoom_out_samples)

    # # run experiment #5 : dataset = zoom in, video length = 16
    # config.ds_type = 'zoom_in'
    # config.sample_len = 16
    # vvn.run(config=config, indices=zoom_in_samples)

    # # run experiment #6 : dataset = zoom out, video length = 16
    # config.ds_type = 'zoom_out'
    # vvn.run(config=config, indices=zoom_out_samples)

