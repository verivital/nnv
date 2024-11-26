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
	
# -------------------------------------------------------------------------
# NNV INSTALLATION
matlab -nodisplay -r "cd ${CURR_DIR}; cd code/nnv/; install; exit()"


# DONE
