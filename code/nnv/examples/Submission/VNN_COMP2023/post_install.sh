cd ~/.matlab/R2022b_licenses

curl --retry 100 --retry-connrefused -L -O 
sleep 60
ls -al

unzip  *.zip*

ls -al

cp -f license.lic /usr/local/matlab/licenses/

rm /usr/local/matlab/licenses/license_info.xml

sudo mkdir /usr/local/MATLAB/R2022b

sudo chmod u+rwx -R /usr/local/MATLAB/R2022b

matlab -nodisplay -r "cd /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2023/; prepare_run; quit"

sudo apt install -y python3-pip 
pip install numpy 

cd /usr/local/matlab/extern/engines/python
python3 -m pip install .

# TEST IF MATLAB ENGINE IS INSTALLED PROPERLY
# START_ENGINE ='import matlab.engine\nimport time\n eng = matlab.engine.start_matlab() \n print(eng) \n eng.prepare_run() \n exit()'
python3 -c "exec('import matlab.engine\nimport time\n eng = matlab.engine.start_matlab() \n print(eng) \n eng.prepare_run() \n exit()')"

# TEST IF WE CAN FIND MATLAB
cd /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2023/
./prepare_instance.sh v1 acasxu acas xu