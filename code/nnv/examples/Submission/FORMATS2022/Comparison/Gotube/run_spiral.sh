#!/bin/bash

# 2 dimensional spiral (linear)
BENCHMARK=spiralL
TIME_STEP=0.01
TIME_HORIZON=10

BATCH_SIZE=1000
timeout 1800 python main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.01 --gamma 0.01 --score --profile
timeout 1800 python main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.05 --gamma 0.01 --score --profile
timeout 1800 python main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.1 --gamma 0.01 --score --profile
