import os
import sys

import vvn.verify as vvn
from vvn.config import Config


if __name__=="__main__":
    # get the results dir
    # root = os.path.dirname(os.getcwd())
    # output_dir = os.path.join(root, 'FORMALISE2025', 'results')
    
    # for SoSym (in docker)
    output_dir = "/tmp"

    # define the starting configurations
    config = Config(
        class_size=11,
        epsilon=[1/255,2/255,3/255],
        ds_type='ucf11',
        sample_len=8,
        ver_algorithm='relax',
        timeout=1800,
        output_dir=output_dir
    )

    #
    # 8 frames
    # 
    # samples = list(range(1, 10))
    # samples = list(range(10, 16))
    samples = list(range(16, 26))
    vvn.run_ucf11(config=config, indices=samples, out_channels=8)
    vvn.run_ucf11(config=config, indices=samples, out_channels=16)
    vvn.run_ucf11(config=config, indices=samples, out_channels=32)
    vvn.run_ucf11(config=config, indices=samples, out_channels=64)    
    # vvn.run_ucf11(config=config, indices=samples, out_channels=128)

    #
    # 16 frames
    #
    config.sample_len = 16
    vvn.run_ucf11(config=config, indices=samples, out_channels=8)
    vvn.run_ucf11(config=config, indices=samples, out_channels=16)
    vvn.run_ucf11(config=config, indices=samples, out_channels=32)
    # vvn.run_ucf11(config=config, indices=samples, out_channels=64)    
    # vvn.run_ucf11(config=config, indices=samples, out_channels=128)


    #
    # 32 frames
    #
    config.sample_len = 32
    vvn.run_ucf11(config=config, indices=samples, out_channels=8)
    vvn.run_ucf11(config=config, indices=samples, out_channels=16)
