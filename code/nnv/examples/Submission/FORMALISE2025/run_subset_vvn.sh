#!/bin/bash
# TABLE II (GET ALL VERIFICATION RESULTS)
python src/run.py subset # MNIST Video
python src/run_gtsrb.py subset # GTSRB
python src/run_stmnist.py subset # STMNIST
python src/analysis/make_table2.py


# FIGURE 8 (COMPARISON OF AVERAGE RUNTIME)
python src/analysis/make_plots.py


# FIGURE 7 (REACHABLE OUTPUT PLOTS)
python src/generate_output_plots.py