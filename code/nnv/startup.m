fprintf('\nAdding dependencies to Matlab path...\n');
tbxmanager restorepath

fprintf('\nAdding NNV to Matlab path...\n');
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
p = genpath(newdir); % generate a path that includes NNV folder and all folders below it
addpath(p);
cd(newdir);
clear;