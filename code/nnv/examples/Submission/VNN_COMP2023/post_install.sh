cd ~/.matlab/R2022b_licenses

curl --retry 100 --retry-connrefused -L -O https://www.dropbox.com/scl/fi/dzu7533mj5lh6ur7dspp2/otherFiles.zip?rlkey=ee03nzds3zog73nw46hmbc510&dl=0
sleep 60
ls -al

unzip  *.zip*

ls -al

cp -f license.lic /usr/local/matlab/licenses/

rm /usr/local/matlab/licenses/license_info.xml

matlab -nodisplay -r "cd /home/ubuntu/toolkit/code/nnv/examples/Submission/VNN_COMP2023/; prepare_run; quit"

sudo apt install -y python3-pip 
pip install numpy 

cd /usr/local/matlab/extern/engines/python
python -m pip install .