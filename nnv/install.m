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

fprintf('\nInstalling MatConvNet....');
mydir  = pwd;
idcs   = strfind(mydir,filesep);
newdir = mydir(1:idcs(end)-1);
cd(newdir);
cd('engine');
cd('matconvnet');
cd('matconvnet-1.0-beta25');
% installing required compiler
% see here: https://developer.microsoft.com/en-us/windows/downloads/windows-10-sdk
% some error may occur: 
% see here: https://www.mathworks.com/matlabcentral/answers/331523-unable-to-find-cl-exe-executing-vl_compilenn
% note: MinGW64 is not supported in matconvnet

% install visual studio 2015 
% here: https://stackoverflow.com/questions/44290672/how-to-download-visual-studio-community-edition-2015-not-2017
% probem may occur after installing windows 10 sdk
% cl.exe problem: see here: https://social.msdn.microsoft.com/Forums/vstudio/en-US/2ebc59de-dc8e-426f-b1ca-885c2709df5f/visual-c-installed-but-cannot-find-clexe?forum=vcgeneral

mex -setup C++;
run matlab/vl_compilenn;

startup_nnv; % adding dependencies and nnv to the path

fprintf('\nInstalling NNV is done, it is ready to use.');
fprintf('\nPlease go to examples or test folders to run case studies and test examples.');
fprintf('\nTHANK YOU FOR TRYING NNV. PLEASE EMAIL trhoangdung@gmail.com FOR ANY CONCERNS OR SUGESSTIONS\n\n');
