tbxmanager_folder = 'tbxmanager';
fprintf('\nUNINSTALLING NNV...');
cd(tbxmanager_folder);
% uninstall mpt toolbox and dependencies
fprintf('\nUninstalling mpt toolbox and dependencies...');
tbxmanager uninstall mpt mptdoc cddmex fourier glpkmex hysdel lcp yalmip sedumi;
fprintf('\nmpt toolbox and dependencies have been removed');

% clean nnv in Matlab Path
fprintf('\nRemoving NNV from Matlab Path');
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
p = genpath(newdir); % generate a path that includes NNV folder and all folders below it
rmpath(p); % remove all folders related to NNV
fprintf('\nNNV has been removed from Matlab search path.');

% delete tbxmanager
fprintf('\nRemoving tbxmanager toolbox...');
cd(newdir);
rmdir tbxmanager s;
fprintf('\ntbxmanager has been removed.');
fprintf('\nUNINSTALLATION IS DONE!. THANK YOU FOR TRYING NNV. PLEASE EMAIL TO trhoangdung@gmail.com for ANY CONCERNS OR SUGGESTIONS');
clear;


