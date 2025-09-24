#!/bin/bash

# VARS
MATLAB_RELEASE=2024b # latest
EXISTING_MATLAB_LOCATION=$(dirname $(dirname $(readlink -f $(which matlab))))
# (assume no products are installed yet)
ADDITIONAL_PRODUCTS="Computer_Vision_Toolbox Control_System_Toolbox Deep_Learning_Toolbox Image_Processing_Toolbox Optimization_Toolbox Parallel_Computing_Toolbox Statistics_and_Machine_Learning_Toolbox Symbolic_Math_Toolbox System_Identification_Toolbox Deep_Learning_Toolbox_Converter_for_ONNX_Model_Format Deep_Learning_Toolbox_Converter_for_TensorFlow_Models"
CURR_DIR=$(pwd)

	
# MATLAB INSTALLATION
wget -q https://www.mathworks.com/mpm/glnxa64/mpm \
    && chmod +x mpm \
    && ./mpm install \
        --destination=${EXISTING_MATLAB_LOCATION} \
        --release=r${MATLAB_RELEASE} \
        --products ${ADDITIONAL_PRODUCTS}	

# -----------------------------------------
# Install python environment for CP
#------------------------------------------
# Define the name of the virtual environment directory
VENV_DIR=".venv"
# Create the virtual environment if it doesn't exist
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating virtual environment: $VENV_DIR"
    python3 -m venv "$VENV_DIR"
else
    echo "Virtual environment $VENV_DIR already exists."
fi

# Activate the virtual environment
echo "Activating virtual environment: $VENV_DIR"
source "$VENV_DIR/bin/activate"

# Upgrade pip within the virtual environment
echo "Upgrading pip..."
pip install --upgrade pip

# Install required packages from a requirements.txt file (optional)
echo "Installing packages from requirement.txt..."
pip install -r requirements.txt



# -------------------------------------------------------------------------
# NNV INSTALLATION
matlab -nodisplay -r "cd ${CURR_DIR}; cd code/nnv/; install; exit()"


# DONE
