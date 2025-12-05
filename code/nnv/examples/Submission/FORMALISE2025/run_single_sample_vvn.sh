#!/bin/bash
# Run a single sample from each dataset variation for quick testing
# This verifies all dataset/algorithm combinations work correctly

OUTPUT_FILE="test_single_all_samples.txt"

# Clear the file if it exists
> $OUTPUT_FILE

echo "Starting single-sample test for all datasets..."
echo "Starting single-sample test for all datasets..." >> $OUTPUT_FILE

#######################
### ZOOM IN RESULTS ###
#######################
echo "Testing Zoom In..."
# Zoom In 4f, relax
python src/vvn/verify.py zoom_in relax 1 4 1800 751 >> $OUTPUT_FILE 2>&1

# Zoom In 8f, relax
python src/vvn/verify.py zoom_in relax 1 8 1800 751 >> $OUTPUT_FILE 2>&1

# Zoom In 16f, relax
python src/vvn/verify.py zoom_in relax 1 16 1800 751 >> $OUTPUT_FILE 2>&1

# Zoom In 4f, approx
python src/vvn/verify.py zoom_in approx 1 4 1800 751 >> $OUTPUT_FILE 2>&1

# Zoom In 8f, approx
python src/vvn/verify.py zoom_in approx 1 8 1800 751 >> $OUTPUT_FILE 2>&1

# Zoom In 16f, approx
python src/vvn/verify.py zoom_in approx 1 16 1800 751 >> $OUTPUT_FILE 2>&1


########################
### ZOOM OUT RESULTS ###
########################
echo "Testing Zoom Out..."
# Zoom Out 4f, relax
python src/vvn/verify.py zoom_out relax 1 4 1800 624 >> $OUTPUT_FILE 2>&1

# Zoom Out 8f, relax
python src/vvn/verify.py zoom_out relax 1 8 1800 624 >> $OUTPUT_FILE 2>&1

# Zoom Out 16f, relax
python src/vvn/verify.py zoom_out relax 1 16 1800 624 >> $OUTPUT_FILE 2>&1

# Zoom Out 4f, approx
python src/vvn/verify.py zoom_out approx 1 4 1800 624 >> $OUTPUT_FILE 2>&1

# Zoom Out 8f, approx
python src/vvn/verify.py zoom_out approx 1 8 1800 624 >> $OUTPUT_FILE 2>&1

# Zoom Out 16f, approx
python src/vvn/verify.py zoom_out approx 1 16 1800 624 >> $OUTPUT_FILE 2>&1


#######################
#### GTSRB RESULTS ####
#######################
echo "Testing GTSRB..."
# GTSRB 4f, relax
python src/vvn/verify.py gtsrb relax 1 4 1800 12597 >> $OUTPUT_FILE 2>&1

# GTSRB 8f, relax
python src/vvn/verify.py gtsrb relax 1 8 1800 12597 >> $OUTPUT_FILE 2>&1

# GTSRB 16f, relax
python src/vvn/verify.py gtsrb relax 1 16 1800 12597 >> $OUTPUT_FILE 2>&1

# GTSRB 4f, approx
python src/vvn/verify.py gtsrb approx 1 4 1800 12597 >> $OUTPUT_FILE 2>&1

# GTSRB 8f, approx
python src/vvn/verify.py gtsrb approx 1 8 1800 12597 >> $OUTPUT_FILE 2>&1

# GTSRB 16f, approx
python src/vvn/verify.py gtsrb approx 1 16 1800 12597 >> $OUTPUT_FILE 2>&1


#######################
### STMNIST RESULTS ###
#######################
echo "Testing ST-MNIST..."
# STMNIST 16f, relax
python src/vvn/verify.py stmnist relax 1 16 1800 311 >> $OUTPUT_FILE 2>&1

# STMNIST 32f, relax
python src/vvn/verify.py stmnist relax 1 32 1800 311 >> $OUTPUT_FILE 2>&1

# STMNIST 64f, relax
python src/vvn/verify.py stmnist relax 1 64 1800 311 >> $OUTPUT_FILE 2>&1

# STMNIST 16f, approx
python src/vvn/verify.py stmnist approx 1 16 1800 311 >> $OUTPUT_FILE 2>&1

# STMNIST 32f, approx
python src/vvn/verify.py stmnist approx 1 32 1800 311 >> $OUTPUT_FILE 2>&1

# STMNIST 64f, approx
python src/vvn/verify.py stmnist approx 1 64 1800 311 >> $OUTPUT_FILE 2>&1


##########################
### KTH ACTIONS RESULTS ###
##########################
echo "Testing KTH Actions..."
# KTH Actions 8f, relax (only relax for KTH)
python src/vvn/verify.py kthactions relax 1 8 1800 0 >> $OUTPUT_FILE 2>&1

# KTH Actions 16f, relax
python src/vvn/verify.py kthactions relax 1 16 1800 0 >> $OUTPUT_FILE 2>&1

# KTH Actions 32f, relax
python src/vvn/verify.py kthactions relax 1 32 1800 0 >> $OUTPUT_FILE 2>&1


#######################
### UCF11 RESULTS ###
#######################
echo "Testing UCF11..."
# UCF11 8f, relax (only relax for UCF11)
python src/vvn/verify.py ucf11 relax 1 8 1800 0 >> $OUTPUT_FILE 2>&1

# UCF11 16f, relax
python src/vvn/verify.py ucf11 relax 1 16 1800 0 >> $OUTPUT_FILE 2>&1

# UCF11 32f, relax
python src/vvn/verify.py ucf11 relax 1 32 1800 0 >> $OUTPUT_FILE 2>&1


echo -e "\n\n**********************************************"
echo "             Single-sample test complete.         "
echo "**********************************************"
echo ""
echo "Results written to: $OUTPUT_FILE"
