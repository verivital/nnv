function prepare_run()

%% Print out some info to prepare for the comp

% remove paths from any prior installation (if any)
% rmpath(genpath('/home/ubuntu/toolkit/code/nnv/')); savepath;

% installing nnv
cd /home/ubuntu/toolkit/code/nnv/;
install;
% save path to NNV
addpath(genpath('/home/ubuntu/toolkit/code/nnv/')); savepath;

%disp("");
%disp("Support package path");
%disp(matlabshared.supportpkg.getSupportPackageRoot);

%disp(" ")
%disp("MATLAB root");
%disp(matlabroot);

% set support package root to import onnx location
%matlabshared.supportpkg.setSupportPackageRoot('/usr/local/MATLAB/R2022b'); 
%addpath(genpath('/usr/local/MATLAB')); 

disp("Support package path");
disp(matlabshared.supportpkg.getSupportPackageRoot);

savepath; 
% quit;

end