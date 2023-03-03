#!/bin/bash

# FPA CTRNN
BENCHMARK=fpaCTRNN
TIME_STEP=0.01
TIME_HORIZON=10
INITIAL_RADIUS=0.01
# new parameters - to be defined
BATCH_SIZE=1000
timeout 300 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score --profile
timeout 300 python3 main.py --time_horizon 2.5 --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score --profile
timeout 300 python3 main.py --time_horizon 0.5 --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score --profile

# 2 dimensional spiral (linear)
BENCHMARK=spiralL
TIME_STEP=0.01
TIME_HORIZON=10

BATCH_SIZE=1000
timeout 300 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.01 --gamma 0.01 --score --profile
timeout 300 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.05 --gamma 0.01 --score --profile
timeout 300 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.1 --gamma 0.01 --score --profile
