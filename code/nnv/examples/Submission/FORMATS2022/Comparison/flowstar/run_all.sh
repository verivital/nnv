#!/bin/bash
mkdir outputs
{ timeout 7200 ./flowstar < Benchmarks/spiral_NL1.model ; } > outputs/spiral_NL1.txt
{ timeout 7200 ./flowstar < Benchmarks/spiral_NL2.model ; } > outputs/spiral_NL2.txt
{ timeout 7200 ./flowstar < Benchmarks/spiral_NL3.model ; } > outputs/spiral_NL3.txt
{ timeout 7200 ./flowstar < Benchmarks/spiral_L1.model ; } > outputs/spiral_L1.txt
{ timeout 7200 ./flowstar < Benchmarks/spiral_L2.model ; } > outputs/spiral_L2.txt
{ timeout 7200 ./flowstar < Benchmarks/spiral_L3.model ; } > outputs/spiral_L3.txt
{ timeout 7200 ./flowstar < Benchmarks/CTRNN_DFP.model ; } > outputs/dfp_time.txt
{ timeout 7200 ./flowstar < Benchmarks/CTRNN_Cartpole.model ; } > outputs/cartpole_time.txt
{ timeout 7200 ./flowstar < Benchmarks/CTRNN_DFP_short.model ; } > outputs/dfp_time_short.txt
{ timeout 7200 ./flowstar < Benchmarks/CTRNN_Cartpole_short.model ; } > outputs/cartpole_time_short.txt
{ timeout 7200 ./flowstar < Benchmarks/CTRNN_DFP_mid.model ; } > outputs/dfp_time_mid.txt
{ timeout 7200 ./flowstar < Benchmarks/CTRNN_Cartpole_mid.model ; } > outputs/cartpole_time_mid.txt

