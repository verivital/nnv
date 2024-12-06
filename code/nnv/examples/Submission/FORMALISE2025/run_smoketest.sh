#!/bin/bash

OUTPUT_FILE="smoketest_outputs.txt"

# Clear the file if it exists
> $OUTPUT_FILE

echo "Starting smoke test..."
echo "Starting smoke test..." >> $OUTPUT_FILE

#######################
### ZOOM IN RESULTS ### 
#######################
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


########################
### ZOOM OUT RESULTS ### 
########################
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


#######################
#### GTSRB RESULTS ####
#######################
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


echo -e "\n\n**********************************************"
echo "              Smoke test complete.         "
echo -e "**********************************************\n\n"
