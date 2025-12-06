#!/bin/bash
# Reproduce a subset of results (first row of each table)
# Useful for quick validation that everything works

echo "=============================================="
echo "Running subset of experiments for validation"
echo "=============================================="

# MNIST Video - first experiment only (Zoom In, 4 frames)
echo ""
echo "Running MNIST Video subset (--subset flag)..."
python src/run.py --subset

# Generate reachable output plots (Figure 7)
echo ""
echo "Generating reachable output plots..."
python src/generate_output_plots.py

echo ""
echo "=============================================="
echo "Subset experiments complete."
echo "=============================================="
