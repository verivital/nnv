NNV AWS installation steps:

1) Set up reference architecture from cloudformation (launch stack), on us-east-1

https://github.com/mathworks-ref-arch/matlab-on-aws/tree/master/releases/R2021a

2) Login via RDP (this step requires GUI to activate from our testing, maybe possible without, but didn't manage to figure out how)

3) Activate Matlab: this will enable terminal usage of matlab command (for no gui mode with matlab -nodisplay)

Terminal (within RDP session as a GUI will get launched) to:
/usr/local/matlab/bin
./activate_matlab.sh
enter matlab account information for a valid user (configured with taylor's account for vnn-comp'21)

Update Matlab if possible (via RDP GUI interface to update to matlab 2021a update 3; not necessarily required)

4) install matlab python engine, details here: https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html

(sudo needed because matlab is installed as user root; if this step is moved, need to change sudo, e.g. if after the chown as this would cause permission problems with some files root and some not)

cd /usr/local/matlab/extern/engines/python
sudo python3 setup.py install

5) install onnx support package (via silent installation from terminal, OR via RDP using typical GUI interaction for support packages)
note: this requires root, however, matlab had an error launching with sudo due to some of the xhost/rdp set up

so, to install support packages, as the user will not be root (but root is required to modify the matlab installation files), need to chown the matlab installation files to ubuntu username/group
other support packages (pytorch, tensorflow, etc.) can be installed if needed

/usr/local$ sudo chown -R ubuntu:ubuntu matlab

Details on installing support packages can be found e.g. here:

https://www.mathworks.com/help/deeplearning/ref/importonnxnetwork.html

Silent mode requires having downloaded the matlab support packages elsewhere, then running a script install file (the support packages are basically just copied to a directory within the matlab installation)

Apparent alternative solutions available in comments from Mathworks support from here:

"Amanjit Dulai
2 Jul 2021

@Taylor Johnson: As a workaround you can use the following command (run in a terminal on the remote machine) to set the ownership such that the "ubuntu" user has correct permissions on the SupportPackage folder:

sudo chown -R ubuntu /usr/local/matlab/SupportPackage/

this should allow you to install support packages via the add-on explorer (running MATLAB without sudo).

Alternatively, you can allow root user access to xserver, by running the following command:

xhost local:root

Then you can use sudo matlab to run MATLAB as root user, which would also enable support packags to be installed via the add-on explorer. It's not usually recommended running GUI applications as sudo (and it sounds like you've already run into issues with this), so the first method is probably preferable."

https://www.mathworks.com/matlabcentral/fileexchange/67296-deep-learning-toolbox-converter-for-onnx-model-format

6) Close matlab in GUI RDP session (as next steps may modify path, etc.)

7) install python dependencies; matlab reference architecture has python 3, but not pip

sudo apt install python3-pip
pip install numpy onnxruntime onnx scipy

numpy
scipy
onnxruntime
onnx

8) run nnv set up (from a terminal)

matlab -nodisplay
cd to nnv directory with install.m
cd work/nnv/code/nnv/
install
savepath
exit

% prior savepath should work without being sudo'd, can check by closing matlab, then starting again and displaying path
% this will show whether all nnv aspects are set up appropriately

9) check the path:
matlab -nodisplay
path

% ensure nnv is on the path (e.g. nnv/code/nnv/engine, etc, should be with prior savepath if permissions allowed modifying matlab path)

10) run basic check

cd ~/work/nnv/code/nnv/examples/Submission/VNN_COMP2021/vnncomp_scripts/

ubuntu@ip-172-31-6-250:~/work/nnv/code/nnv/examples/Submission/VNN_COMP2021/vnncomp_scripts$ chmod 755 demo_test_instances.sh
ubuntu@ip-172-31-6-250:~/work/nnv/code/nnv/examples/Submission/VNN_COMP2021/vnncomp_scripts$ ./demo_test_instances.sh

*) recommend updating matlab to 2021a release 3 (done within GUI)


% other minor note: the nnv clone used by vnn-comp scripts may not work right, as nnv has submodules, so needed to be cloned recursively (will check), as things like cora not needed (vnn-comp wasn't cloned recursively, but should be okay, as main submodule dependency is cora, which isn't used here)
