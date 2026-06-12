cd ~/.matlab/R2024b_licenses

curl --retry 100 --retry-connrefused -L -O https://www.dropbox.com/scl/fi/94gyrz1sdpiuvd7e8nc0r/LicenseTL.zip?rlkey=cbbrrkley11nnxd5r9mgrnkg1&st=4bx85clk&dl=0
sleep 60
ls -al

unzip  *.zip*

ls -al

cp -f license.lic /usr/local/matlab/licenses/

rm *.zip*

rm /usr/local/matlab/licenses/license_info.xml

cd /usr/local/matlab/extern/engines/python
python3 -m pip install .

# ADD STEPS TO INSTALL GUROBI

#cp -f gurobi.lic ~/gurobi1102/

#echo 'export GUROBI_HOME="~/gurobi1102/linux64"' >> ~/.bashrc
#echo 'export GRB_LICENSE_FILE="~/gurobi1102/gurobi.lic"' >> ~/.bashrc
#echo 'export PATH="${PATH}:${GUROBI_HOME}/bin"' >> ~/.bashrc
#echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"' >> ~/.bashrc

# Ensure installation is correct
matlab -nodisplay -r "cd /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2025/; prepare_run; quit"

sudo apt install -y python3-pip
pip install torch
pip install numpy
pip install scipy
# onnxruntime backs BOTH the SAT-witness replay gate (validate_witness_onnx.m ->
# onnx_replay_check.py) and the collins_aerospace falsification path
# (collins_falsify.py). Without it those degrade safely to unknown -- but we'd
# leave every collins point and the -150 witness guard on the table.
pip install onnx onnxruntime

# For next year, let's fix gpu drivers to ensure no potential errors there...
# Enable GPU persistence mode (prevents driver unloading)
sudo nvidia-smi -pm 1

# Lock the kernenl verison and GPU drivers.
sudo apt-mark hold linux-image-generic linux-headers-generic nvidia-driver-535
sudo systemctl disable unattended-upgrades