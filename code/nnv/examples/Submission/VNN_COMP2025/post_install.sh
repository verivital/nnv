#cd ~/.matlab/R2024b_licenses

# curl --retry 100 --retry-connrefused -L -O 
# sleep 60
# ls -al

#unzip  *.zip*

#ls -al

# cp -f license.lic /usr/local/matlab/licenses/

# rm *.zip*

rm /usr/local/matlab/licenses/license_info.xml

export MLM_LICENSE_FILE=

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
pip install numpy 

# TEST IF MATLAB ENGINE IS INSTALLED PROPERLY
# START_ENGINE ='import matlab.engine\nimport time\neng = matlab.engine.start_matlab() \nprint(eng) \neng.prepare_run() \nexit()'
# python3 -c "exec('import matlab.engine\nimport time\neng = matlab.engine.start_matlab() \nprint(eng) \neng.prepare_run(nargout = 0, background=True) \nexit()')"

# TEST IF WE CAN FIND MATLAB
# cd /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2025/
# ./prepare_instance.sh v1 acasxu acas xu

