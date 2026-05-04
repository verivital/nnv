#!/bin/bash
# Run all verification experiments for SoSym paper
# Tables 3, 4, 5

echo "=============================================="
echo "Running all verification experiments"
echo "=============================================="

# Tables 3 & 4: MNIST Video (Zoom In/Out) - relax and approx
echo ""
echo "Running MNIST Video experiments..."
python src/run.py

# Tables 3 & 4: GTSRB - relax and approx
echo ""
echo "Running GTSRB experiments..."
python src/run_gtsrb.py

# Tables 3 & 4: ST-MNIST - relax and approx
echo ""
echo "Running ST-MNIST experiments..."
python src/run_stmnist.py

# Table 3: KTH Actions - relax only
echo ""
echo "Running KTH Actions experiments..."
python src/run_kthactions.py

# Table 5: UCF11 - relax only
echo ""
echo "Running UCF11 experiments..."
python src/run_ucf11.py

# Generate analysis outputs
# python src/analysis/make_table2.py
# python src/analysis/make_plots.py
# python src/generate_output_plots.py

echo ""
echo "=============================================="
echo "All experiments complete."
echo "Results saved to /tmp/results/"
echo "=============================================="
