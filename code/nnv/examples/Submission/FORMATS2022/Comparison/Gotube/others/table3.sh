#!/bin/bash

source venv/bin/activate


TIME_STEP=0.5
TIME_HORIZON=2
# new parameters - to be defined
BATCH_SIZE=2000
python main.py --time_horizon $TIME_HORIZON --mu 1.6 --benchmark CTRNNosc --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.001 --gamma 0.2 --score
python main.py --time_horizon $TIME_HORIZON --mu 1.6 --benchmark CTRNNosc --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.005 --gamma 0.2 --score
python main.py --time_horizon $TIME_HORIZON --mu 1.6 --benchmark CTRNNosc --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.01 --gamma 0.2 --score


TIME_STEP=0.5
TIME_HORIZON=10
# new parameters - to be defined
BATCH_SIZE=5000
python main.py --time_horizon $TIME_HORIZON --benchmark pendulumCTRNN --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.01 --gamma 0.05 --score
python main.py --time_horizon $TIME_HORIZON --benchmark pendulumCTRNN --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.05 --gamma 0.05 --score
python main.py --time_horizon $TIME_HORIZON --benchmark pendulumCTRNN --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.1 --gamma 0.05 --score

TIME_STEP=0.5
TIME_HORIZON=10
# new parameters - to be defined
BATCH_SIZE=5000
python main.py --time_horizon $TIME_HORIZON --benchmark ldsCTRNN --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.01 --gamma 0.05 --score
python main.py --time_horizon $TIME_HORIZON --benchmark ldsCTRNN --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.05 --gamma 0.05 --score
python main.py --time_horizon $TIME_HORIZON --benchmark ldsCTRNN --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.1 --gamma 0.05 --score
