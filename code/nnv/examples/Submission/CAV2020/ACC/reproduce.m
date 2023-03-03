path_reproduce = pwd();
%path_nnv_root = ['..', filesep, '..', filesep, '..', filesep]; % in main folder

path_nnv_root = ['..', filesep, '..', filesep, '..', filesep, '..', filesep]; % in acc subfolder

cd(path_nnv_root);
% run installation if not on codeocean / not already set up
try
    is_codeocean();
catch
    install;
end

% check if CORA is set up, if submodules are not pulled recursively, it won't be
% so ensure you've pulled with: 
% git clone --recursive https://github.com/verivital/nnv.git
try
    coraroot();
catch
    'ERROR: CORA not detected, ensure submodules have been pulled, CORA is required for nonlinear plant analysis'
    return;
end

cd(path_reproduce);

mkdir([path_results(), filesep, 'ACC']);

% set number of cores based on system availability, up to a maximum of 16
config_parallelism(16);

'Figure 4: Generating ACC with linear plant model reachable sets'
plot_linear_ACC_reachSets;

'Table 3 (linear plant): Generating ACC results for linear plant model'
verify_linear_ACC;

'Table 3 (nonlinear plant): Generating ACC results for nonlinear plant model'
verify_nonlinear_ACC;
