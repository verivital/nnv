function prepare_run()

%% Print out some info to prepare for the comp

% remove paths from any prior installation (if any)
% rmpath(genpath('/home/ubuntu/toolkit/code/nnv/')); savepath;

% install gurobi
%cd ~/gurobi1102/linux64/matlab;
%gurobi_setup;

% installing nnv
cd /home/ubuntu/toolkit/code/nnv/;
install;
% save path to NNV
addpath(genpath('/home/ubuntu/toolkit/code/nnv/')); savepath;

disp("MATLAB root");
disp(matlabroot);

disp(" ")

disp("Support package path");
disp(matlabshared.supportpkg.getSupportPackageRoot);

% matlabshared.supportpkg.getInstalled

% savepath; 
% quit;

end