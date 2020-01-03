fprintf('\nINSTALLING NNV....');
fprintf('\nIntalling tbxmanager (requires Matlab R2009a or later) ...');
tbxmanager_folder = 'tbxmanager';

root_folder = pwd();

list = dir;
if ~isfolder(tbxmanager_folder)
    mkdir(tbxmanager_folder);
end

% install mpt toobox and other dependencies
cd(tbxmanager_folder);
urlwrite('https://raw.githubusercontent.com/verivital/tbxmanager/master/tbxmanager.m', 'tbxmanager.m');
%urlwrite('http://www.tbxmanager.com/tbxmanager.m', 'tbxmanager.m');
tbxmanager
savepath
fprintf('\nInstalling tbxmanager toolbox is done!');

cd(root_folder);

fprintf('\nIntalling mpt toolbox and other dependencies...\n');
tbxmanager install mpt mptdoc;
tbxmanager install lcp hysdel cddmex clpmex glpkmex fourier sedumi;
tbxmanager install yalmip; % todo: error due to license, need to force acceptance
fprintf('\nInstalling dependencies is done!');
adjust_glpk;

startup_nnv; % adding dependencies and nnv to the path

fprintf('\nInstalling NNV is done, it is ready to use.');
fprintf('\nPlease go to examples or test folders to run case studies and test examples.');
fprintf('\nTHANK YOU FOR TRYING NNV. PLEASE EMAIL trhoangdung@gmail.com FOR ANY CONCERNS OR SUGESSTIONS\n\n');

function adjust_glpk()
    fid = fopen('../tbxmanager/toolboxes/glpkmex/1.0/glnxa64/glpkmex_1_0_glnxa64/glpk.m');
    cac = textscan( fid, '%s', 'Delimiter','\n','whitespace', '');
    fclose(fid)
    fid = fopen('../tbxmanager/toolboxes/glpkmex/1.0/glnxa64/glpkmex_1_0_glnxa64/glpk.m', 'w');
    change_here = 372;
    for jj = 1 : change_here-1
        fprintf(fid, '%s\n', cac{1}{jj});
    end
    fprintf(fid, '%s\n', '%clear glpkcc;');
    for jj = change_here+1: length(cac{1})
        fprintf(fid, '%s\n', cac{1}{jj});
    end
    fclose(fid);
end
