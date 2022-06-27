#!/bin/bash

source venv/bin/activate

BENCHMARK=cartpoleLTC_RK
#TIME_STEP=1e-6
TIME_STEP=0.01
TIME_HORIZON_SHORT=0.35
TIME_HORIZON_LONG=10
INITIAL_RADIUS=1e-4
# new parameters - to be defined
python main.py --time_horizon $TIME_HORIZON_SHORT --benchmark $BENCHMARK --batch_size 10000 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.05 --score
python main.py --time_horizon $TIME_HORIZON_LONG --benchmark $BENCHMARK --batch_size 10000 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.05 --score


BENCHMARK=cartpoleCTRNN
TIME_STEP=0.1
TIME_HORIZON_SHORT=1
TIME_HORIZON_LONG=10
INITIAL_RADIUS=1e-4
# new parameters - to be defined
python main.py --time_horizon $TIME_HORIZON_SHORT --benchmark $BENCHMARK --batch_size 10000 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.05 --score
python main.py --time_horizon $TIME_HORIZON_LONG --benchmark $BENCHMARK --batch_size 10000 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.05 --score
