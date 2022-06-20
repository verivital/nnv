#!/bin/bash
mkdir short_outputs
{ timeout 300 ./flowstar < Benchmarks/spiral_NL1.model ; } > short_outputs/spiral_NL1.txt
{ timeout 300 ./flowstar < Benchmarks/spiral_NL2.model ; } > short_outputs/spiral_NL2.txt
{ timeout 300 ./flowstar < Benchmarks/spiral_NL3.model ; } > short_outputs/spiral_NL3.txt
{ timeout 300 ./flowstar < Benchmarks/spiral_L1.model ; } > short_outputs/spiral_L1.txt
{ timeout 300 ./flowstar < Benchmarks/spiral_L2.model ; } > short_outputs/spiral_L2.txt
{ timeout 300 ./flowstar < Benchmarks/spiral_L3.model ; } > short_outputs/spiral_L3.txt
{ timeout 300 ./flowstar < Benchmarks/CTRNN_FPA.model ; } > short_outputs/fpa_time.txt
{ timeout 300 ./flowstar < Benchmarks/CTRNN_Cartpole.model ; } > short_outputs/cartpole_time.txt
{ timeout 300 ./flowstar < Benchmarks/CTRNN_FPA_short.model ; } > short_outputs/fpa_time_short.txt
{ timeout 300 ./flowstar < Benchmarks/CTRNN_Cartpole_short.model ; } > short_outputs/cartpole_time_short.txt
{ timeout 300 ./flowstar < Benchmarks/CTRNN_FPA_mid.model ; } > short_outputs/fpa_time_mid.txt
{ timeout 300 ./flowstar < Benchmarks/CTRNN_Cartpole_mid.model ; } > short_outputs/cartpole_time_mid.txt

