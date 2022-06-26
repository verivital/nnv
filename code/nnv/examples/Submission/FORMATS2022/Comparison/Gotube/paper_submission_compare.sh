#!/bin/bash

# #Damped Forced Pendulum CTRNN
BENCHMARK=fpaCTRNN
TIME_STEP=0.01
TIME_HORIZON=10
INITIAL_RADIUS=0.01
# new parameters - to be defined
BATCH_SIZE=1000
timeout 7200 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score --profile
timeout 7200 python3 main.py --time_horizon 2.5 --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score --profile
timeout 7200 python3 main.py --time_horizon 0.5 --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score --profile

# 2 dimensional spiral (linear)
BENCHMARK=spiralL
TIME_STEP=0.01
TIME_HORIZON=10

BATCH_SIZE=1000
timeout 7200 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.01 --gamma 0.01 --score --profile
timeout 7200 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.05 --gamma 0.01 --score --profile
timeout 7200 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.1 --gamma 0.01 --score --profile

# 2 dimensional spiral (nonlinear)
BENCHMARK=spiralNL
TIME_STEP=0.01
TIME_HORIZON=10
BATCH_SIZE=1000
timeout 7200 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.01 --gamma 0.01 --score --profile
timeout 7200 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.05 --gamma 0.01 --score --profile
timeout 7200 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius 0.1 --gamma 0.01 --score --profile


# #Cartpole CTRNN
BENCHMARK=cartpoleCTRNN
TIME_STEP=0.01
TIME_HORIZON=2
INITIAL_RADIUS=1e-4
# new parameters - to be defined
BATCH_SIZE=10000
timeout 7200 python3 main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score --profile
timeout 7200 python3 main.py --time_horizon 1.0 --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score --profile
timeout 7200 python3 main.py --time_horizon 0.1 --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score --profile

