fprintf('\nINSTALLING NNV....');
fprintf('\nIntalling tbxmanager (requires Matlab R2009a or later) ...');
tbxmanager_folder = 'tbxmanager';
list = dir;
if ~isfolder(tbxmanager_folder)
    mkdir(tbxmanager_folder);
end

% install mpt toobox and other dependencies
cd(tbxmanager_folder);
urlwrite('http://www.tbxmanager.com/tbxmanager.m', 'tbxmanager.m');
tbxmanager
savepath
fprintf('\nInstalling tbxmanager toolbox is done!');

fprintf('\nIntalling mpt toolbox and other dependencies...\n');
tbxmanager install mpt mptdoc;
tbxmanager install lcp hysdel cddmex clpmex glpkmex fourier sedumi yalmip;
fprintf('\nInstalling dependencies is done!');

mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
cd(newdir);
startup_nnv; % adding dependencies and nnv to the path

fprintf('\nInstalling NNV is done, it is ready to use.');
fprintf('\nPlease go to examples or test folders to run case studies and test examples.');
fprintf('\nTHANK YOU FOR TRYING NNV. PLEASE EMAIL trhoangdung@gmail.com FOR ANY CONCERNS OR SUGESSTIONS\n\n');
