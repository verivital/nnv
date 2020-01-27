path_reproduce = pwd();
%path_nnv_root = ['..', filesep, '..', filesep, '..', filesep]; % in main folder

path_nnv_root = ['..', filesep, '..', filesep, '..', filesep, '..', filesep]; % in acc subfolder

cd(path_nnv_root);
install;

cd(path_reproduce);

% set number of cores to use in parallel computations
% if system executed on has fewer, set lower
physical_cores = feature('numcores'); % physical cores

% use at most the number of logical cores available
core_info = evalc('feature(''numcores'')'); % string with all core info, including logical cores
core_info_logical = regexpi(core_info, '(?<num>\d+)\s+logical\s+cores', 'names');
available_cores = physical_cores;
for i = 1 : length(core_info_logical)
    available_cores = max(available_cores, str2num(core_info_logical(i).num));
end

numCores = min(available_cores, 4); % use at most 4 cores if available, else fewer

plot_linear_ACC_reachSets;

verify_linear_ACC;

verify_nonlinear_ACC;
