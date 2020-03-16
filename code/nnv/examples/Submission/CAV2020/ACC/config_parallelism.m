function [out] = config_parallelism(max_cores)
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

    numCores = min(available_cores, max_cores); % use at most max_cores cores if available, else fewer

    out = numCores;
end
