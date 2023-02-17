#!/bin/bash

python main.py --time_horizon 2 --mu 2.0 --benchmark CTRNNosc --batch_size 256 --time_step 0.5 --radius 0.001 --gamma 0.5 --score
python main.py --time_horizon 2 --mu 1.8 --benchmark CTRNNosc --batch_size 256 --time_step 0.5 --radius 0.001 --gamma 0.5 --score
python main.py --time_horizon 2 --mu 1.7 --benchmark CTRNNosc --batch_size 256 --time_step 0.5 --radius 0.001 --gamma 0.5 --score
python main.py --time_horizon 2 --mu 1.6 --benchmark CTRNNosc --batch_size 256 --time_step 0.5 --radius 0.001 --gamma 0.5 --score
python main.py --time_horizon 2 --mu 1.5 --benchmark CTRNNosc --batch_size 1024 --time_step 0.5 --radius 0.001 --gamma 0.5 --score

python main.py --time_horizon 10 --mu 2.0 --benchmark pendulumCTRNN --batch_size 256 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.8 --benchmark pendulumCTRNN --batch_size 256 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.7 --benchmark pendulumCTRNN --batch_size 256 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.6 --benchmark pendulumCTRNN --batch_size 256 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.5 --benchmark pendulumCTRNN --batch_size 1024 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.3 --benchmark pendulumCTRNN --batch_size 1024 --time_step 0.5 --radius 0.01 --gamma 0.5 --score

python main.py --time_horizon 10 --mu 2.0 --benchmark ldsCTRNN --batch_size 256 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.8 --benchmark ldsCTRNN --batch_size 256 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.7 --benchmark ldsCTRNN --batch_size 256 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.6 --benchmark ldsCTRNN --batch_size 256 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.5 --benchmark ldsCTRNN --batch_size 1024 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.3 --benchmark ldsCTRNN --batch_size 1024 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
python main.py --time_horizon 10 --mu 1.1 --benchmark ldsCTRNN --batch_size 1024 --time_step 0.5 --radius 0.01 --gamma 0.5 --score
