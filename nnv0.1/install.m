fprintf('\nINSTALLING NNV....')
fprintf('\nIntalling tbxmanager (requires Matlab R2009a or later) ...');
tbxmanager_folder = 'tbxmanager';
list = dir;
if ~exist(tbxmanager_folder)
    mkdir(tbxmanager_folder);
end

% install mpt toobox and other dependencies
cd(tbxmanager_folder);
urlwrite('http://www.tbxmanager.com/tbxmanager.m', 'tbxmanager.m');
tbxmanager
savepath
tbxmanager restorepath
fprintf('\nInstalling tbxmanager toolbox is done!');

fprintf('\nIntalling mpt toolbox and other dependencies...');
tbxmanager install mpt mptdoc;
tbxmanager install lcp hysdel cddmex clpmex glpkmex fourier sedumi yalmip;
fprintf('\nInstalling dependencies is done!');

% adding nnv to the path
fprintf('\nAdding NNV to Matlab path...');
mydir  = pwd;
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end)-1);
p = genpath(newdir); % generate a path that includes NNV folder and all folders below it
addpath(p);
fprintf('\nInstalling NNV is done, it is ready to use.');
fprintf('\nPlease go to examples or test folders to run case studies and test examples.');
fprintf('\nTHANK YOU FOR TRYING NNV. PLEASE EMAIL trhoangdung@gmail.com FOR ANY CONCERNS OR SUGESSTIONS');
cd(newdir);
clear;
