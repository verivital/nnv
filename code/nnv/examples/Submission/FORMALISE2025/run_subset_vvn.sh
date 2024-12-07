#!/bin/bash

# Reproduce the first row of table 2.

# TABLE II (GET ALL VERIFICATION RESULTS)
python src/run.py --subset # MNIST Video

# FIGURE 7 (REACHABLE OUTPUT PLOTS)
python src/generate_output_plots.py