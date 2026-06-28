mkdir -p ~/.matlab/R2026a_licenses
cd ~/.matlab/R2026a_licenses

# MAC-locked MATLAB R2026a license (HOSTID 02e1e896fadb == ENI eni-0b11771dfe21b94ee).
# Verified 2026-06-28: FLEXlm passcode file, 116 products incl. Deep Learning (Neural_Network_Toolbox),
# Parallel Computing (Distrib_Computing_Toolbox, for GPU/parfor), Optimization (linprog), GPU_Coder;
# expires 30-may-2027. Fetched at runtime; the .lic is MAC-locked so the URL is unusable off this ENI.
# NOTE: URL is QUOTED (the prior year's unquoted &-URL backgrounded curl + dropped query params).
curl --retry 100 --retry-connrefused -L -o license.lic "https://www.dropbox.com/scl/fi/w5jgddmf3qm5znjw67ajm/matlab-license-vnncomp2026-nnv.lic?rlkey=z3wnimbad4ykjyde95yhq7ik3&st=2oc64st3&dl=1"
sleep 5
ls -al

# sudo: /usr/local/matlab/licenses is root-owned and post_install runs non-root
# (run_post_installation_script_as_root=False). Without sudo the copy fails -> MATLAB never
# licenses -> every benchmark unknown. Passwordless sudo is available (install_tool.sh uses it).
sudo cp -f license.lic /usr/local/matlab/licenses/

sudo rm -f /usr/local/matlab/licenses/license_info.xml

cd /usr/local/matlab/extern/engines/python
python3 -m pip install .

# ADD STEPS TO INSTALL GUROBI

#cp -f gurobi.lic ~/gurobi1102/

#echo 'export GUROBI_HOME="~/gurobi1102/linux64"' >> ~/.bashrc
#echo 'export GRB_LICENSE_FILE="~/gurobi1102/gurobi.lic"' >> ~/.bashrc
#echo 'export PATH="${PATH}:${GUROBI_HOME}/bin"' >> ~/.bashrc
#echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"' >> ~/.bashrc

# Ensure installation is correct
matlab -nodisplay -r "cd /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2026/; prepare_run; quit"

sudo apt install -y python3-pip
# torch backs engine/nn/Prob_reach (cp-star / probabilistic-reach paths) -- keep it.
python3 -m pip install torch
python3 -m pip install numpy
python3 -m pip install scipy
# onnxruntime backs BOTH the SAT-witness replay gate (validate_witness_onnx.m ->
# onnx_replay_check.py) and the collins_aerospace falsification path
# (collins_falsify.py). Without it those degrade safely to unknown -- but we'd
# leave every collins point and the -150 witness guard on the table.
# pinned to the versions validated on the 2026-06-12 AWS dry run (m5.16xlarge, R2026a)
python3 -m pip install onnx==1.20.0 onnxruntime==1.23.1

# For next year, let's fix gpu drivers to ensure no potential errors there...
# Enable GPU persistence mode (prevents driver unloading)
sudo nvidia-smi -pm 1

# Lock the kernel version and GPU drivers (driver 570 -- the version on our tested g5 AMI; NNV GPU
# reach needs >=570). Pair with the form's "restart after post-install" so the driver is reloaded.
sudo apt-mark hold linux-image-generic linux-headers-generic nvidia-driver-570
sudo systemctl disable unattended-upgrades