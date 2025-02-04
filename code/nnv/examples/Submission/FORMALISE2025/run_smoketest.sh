#!/bin/bash

echo "Starting smoke test..."

# Zoom In 4f, relax
python src/vvn/verify.py zoom_in relax 1 4 1800 751

echo -e "\n\n**********************************************"
echo "              Smoke test complete.         "
echo -e "**********************************************\n\n"
