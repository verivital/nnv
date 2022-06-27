#!/bin/bash

source venv/bin/activate

#Brusselator
# BENCHMARK=bruss
MU=1.3
TIME_STEP=0.01
TIME_HORIZON=9
INITIAL_RADIUS=0.01
# new parameters - to be defined
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 5 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.5 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 5 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.1 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 5 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score

# #Van der Pol
BENCHMARK=vdp
TIME_STEP=0.01
TIME_HORIZON=40
INITIAL_RADIUS=0.01
# new parameters - to be defined
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 5 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.5 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 5 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.1 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 5 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score

# #Robotarm
BENCHMARK=robot
TIME_STEP=0.01
TIME_HORIZON=40
INITIAL_RADIUS=0.005
# new parameters - to be defined
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 5 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.5 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 5 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.1 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 5 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score

#Dubins car
BENCHMARK=dubins
TIME_STEP=0.1
TIME_HORIZON=60
INITIAL_RADIUS=0.01
# new parameters - to be defined
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 10 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.5 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 10 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.1 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 10 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score

# #MS cardiac cell
BENCHMARK=ms
TIME_STEP=0.01
TIME_HORIZON=10
INITIAL_RADIUS=1e-4
# new parameters - to be defined
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 10 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.5 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 10 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.1 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 10 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score

# #Controlled Cartpole
BENCHMARK=cartpole
TIME_STEP=0.01
TIME_HORIZON=10
INITIAL_RADIUS=1e-4
# new parameters - to be defined
BATCH_SIZE=1000
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 20 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.5 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 20 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.1 --score
python main.py --time_horizon $TIME_HORIZON --mu $MU --benchmark $BENCHMARK --batch_size 20 --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score


# #Cartpole CTRNN
BENCHMARK=cartpoleCTRNN
TIME_STEP=0.02
TIME_HORIZON=1
INITIAL_RADIUS=1e-4
# new parameters - to be defined
BATCH_SIZE=10000
python main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.5 --score
python main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.1 --score
python main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score

# #Cartpole LTC
BENCHMARK=cartpoleLTC_RK
#TIME_STEP=1e-6
TIME_STEP=0.01
TIME_HORIZON=0.35
INITIAL_RADIUS=1e-4
# new parameters - to be defined
BATCH_SIZE=10000
python main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.5 --score
python main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.1 --score
python main.py --time_horizon $TIME_HORIZON --benchmark $BENCHMARK --batch_size $BATCH_SIZE --time_step $TIME_STEP --radius $INITIAL_RADIUS --gamma 0.01 --score

