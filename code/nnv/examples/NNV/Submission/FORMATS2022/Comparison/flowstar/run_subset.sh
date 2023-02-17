#!/bin/bash
mkdir -p /results/logs/flowresults
{ timeout 300 ./flowstar < Benchmarks/spiral_L1.model ; } > /results/logs/flowresults/spiral_L1.txt
{ timeout 300 ./flowstar < Benchmarks/spiral_L2.model ; } > /results/logs/flowresults/spiral_L2.txt
{ timeout 300 ./flowstar < Benchmarks/spiral_L3.model ; } > /results/logs/flowresults/spiral_L3.txt
{ timeout 300 ./flowstar < Benchmarks/CTRNN_FPA.model ; } > /results/logs/flowresults/fpa_time.txt
{ timeout 300 ./flowstar < Benchmarks/CTRNN_FPA_short.model ; } > /results/logs/flowresults/fpa_time_short.txt
{ timeout 300 ./flowstar < Benchmarks/CTRNN_FPA_mid.model ; } > /results/logs/flowresults/fpa_time_mid.txt
mv /flowstar/flowstar-2.1.0/outputs /results/logs/flowresults/
