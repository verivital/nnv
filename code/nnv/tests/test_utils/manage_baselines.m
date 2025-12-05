function results = manage_baselines(action, varargin)
    % MANAGE_BASELINES - Manage test baseline files for regression detection
    %
    % Usage:
    %   manage_baselines('save')      - Save current test data as new baselines
    %   manage_baselines('compare')   - Compare current data against baselines
    %   manage_baselines('list')      - List all baseline files
    %   manage_baselines('status')    - Show baseline status (missing/outdated)
    %   manage_baselines('clean')     - Remove all current test data (not baselines)
    %
    % Options:
    %   manage_baselines('save', 'subdir', 'set')  - Only save baselines in set/
    %   manage_baselines('compare', 'verbose', true) - Show detailed comparison
    %
    % Examples:
    %   % After running tests, save current outputs as baselines
    %   manage_baselines('save');
    %
    %   % In CI/CD, compare test outputs against baselines
    %   results = manage_baselines('compare');
    %   if ~isempty(results.regressions)
    %       error('Regression detected!');
    %   end

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'subdir', '', @ischar);
    addParameter(p, 'verbose', false, @islogical);
    addParameter(p, 'tolerance', 1e-6, @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;

    % Get directories from config
    config = get_test_config();
    data_dir = config.data_dir;
    baselines_dir = config.baselines_dir;

    % Apply subdir filter if specified
    if ~isempty(opts.subdir)
        data_dir = fullfile(data_dir, opts.subdir);
        baselines_dir = fullfile(baselines_dir, opts.subdir);
    end

    results = struct();
    results.action = action;
    results.timestamp = datetime('now');

    switch lower(action)
        case 'save'
            results = save_baselines(data_dir, baselines_dir, opts);

        case 'compare'
            results = compare_baselines(data_dir, baselines_dir, opts);

        case 'list'
            results = list_baselines(baselines_dir);

        case 'status'
            results = check_status(data_dir, baselines_dir);

        case 'clean'
            results = clean_data(data_dir);

        otherwise
            error('Unknown action: %s. Use save, compare, list, status, or clean.', action);
    end
end

function results = save_baselines(data_dir, baselines_dir, opts)
    % Save current test data as new baselines
    results = struct();
    results.saved = {};
    results.errors = {};

    if ~exist(data_dir, 'dir')
        warning('Data directory does not exist: %s', data_dir);
        return;
    end

    % Create baselines directory if needed
    if ~exist(baselines_dir, 'dir')
        mkdir(baselines_dir);
    end

    % Find all .mat files in data directory (recursively)
    mat_files = dir(fullfile(data_dir, '**', '*.mat'));

    for i = 1:length(mat_files)
        src_file = fullfile(mat_files(i).folder, mat_files(i).name);

        % Determine relative path from data_dir
        rel_path = strrep(mat_files(i).folder, data_dir, '');
        if ~isempty(rel_path) && (rel_path(1) == filesep)
            rel_path = rel_path(2:end);
        end

        % Create destination directory
        if isempty(rel_path)
            dest_dir = baselines_dir;
        else
            dest_dir = fullfile(baselines_dir, rel_path);
        end
        if ~exist(dest_dir, 'dir')
            mkdir(dest_dir);
        end

        % Copy file
        dest_file = fullfile(dest_dir, mat_files(i).name);
        try
            copyfile(src_file, dest_file);
            results.saved{end+1} = mat_files(i).name;
            if opts.verbose
                fprintf('  Saved baseline: %s\n', mat_files(i).name);
            end
        catch ME
            results.errors{end+1} = struct('file', mat_files(i).name, 'error', ME.message);
        end
    end

    fprintf('Saved %d baselines to %s\n', length(results.saved), baselines_dir);
end

function results = compare_baselines(data_dir, baselines_dir, opts)
    % Compare current test data against baselines
    results = struct();
    results.matches = {};
    results.regressions = {};
    results.missing_baselines = {};
    results.missing_data = {};

    if ~exist(baselines_dir, 'dir')
        warning('Baselines directory does not exist: %s', baselines_dir);
        warning('Run manage_baselines(''save'') first to create baselines.');
        return;
    end

    % Find all baseline .mat files
    baseline_files = dir(fullfile(baselines_dir, '**', '*.mat'));

    for i = 1:length(baseline_files)
        baseline_file = fullfile(baseline_files(i).folder, baseline_files(i).name);

        % Determine relative path
        rel_path = strrep(baseline_files(i).folder, baselines_dir, '');
        if ~isempty(rel_path) && (rel_path(1) == filesep)
            rel_path = rel_path(2:end);
        end

        % Find corresponding data file
        if isempty(rel_path)
            data_file = fullfile(data_dir, baseline_files(i).name);
        else
            data_file = fullfile(data_dir, rel_path, baseline_files(i).name);
        end

        if ~exist(data_file, 'file')
            results.missing_data{end+1} = baseline_files(i).name;
            if opts.verbose
                fprintf('  Missing data: %s\n', baseline_files(i).name);
            end
            continue;
        end

        % Compare the files
        [match, diff] = compare_mat_files(data_file, baseline_file, opts.tolerance);

        if match
            results.matches{end+1} = baseline_files(i).name;
            if opts.verbose
                fprintf('  Match: %s\n', baseline_files(i).name);
            end
        else
            results.regressions{end+1} = struct(...
                'file', baseline_files(i).name, ...
                'diff', diff);
            fprintf('  REGRESSION: %s\n', baseline_files(i).name);
            if opts.verbose && ~isempty(diff)
                disp(diff);
            end
        end
    end

    % Check for new data files without baselines
    data_files = dir(fullfile(data_dir, '**', '*.mat'));
    for i = 1:length(data_files)
        rel_path = strrep(data_files(i).folder, data_dir, '');
        if ~isempty(rel_path) && (rel_path(1) == filesep)
            rel_path = rel_path(2:end);
        end

        if isempty(rel_path)
            baseline_file = fullfile(baselines_dir, data_files(i).name);
        else
            baseline_file = fullfile(baselines_dir, rel_path, data_files(i).name);
        end

        if ~exist(baseline_file, 'file')
            results.missing_baselines{end+1} = data_files(i).name;
        end
    end

    % Print summary
    fprintf('\nBaseline Comparison Summary:\n');
    fprintf('  Matches: %d\n', length(results.matches));
    fprintf('  Regressions: %d\n', length(results.regressions));
    fprintf('  Missing baselines: %d\n', length(results.missing_baselines));
    fprintf('  Missing data: %d\n', length(results.missing_data));

    if ~isempty(results.regressions)
        fprintf('\n*** REGRESSIONS DETECTED ***\n');
        for i = 1:length(results.regressions)
            fprintf('  - %s\n', results.regressions{i}.file);
        end
    end
end

function [match, diff] = compare_mat_files(data_file, baseline_file, tolerance)
    % Compare two .mat files field by field
    match = true;
    diff = struct();

    try
        data = load(data_file);
        baseline = load(baseline_file);
    catch ME
        match = false;
        diff.error = ME.message;
        return;
    end

    % Get all field names
    data_fields = fieldnames(data);
    baseline_fields = fieldnames(baseline);

    % Check for missing fields
    missing_in_data = setdiff(baseline_fields, data_fields);
    missing_in_baseline = setdiff(data_fields, baseline_fields);

    if ~isempty(missing_in_data) || ~isempty(missing_in_baseline)
        match = false;
        diff.missing_in_data = missing_in_data;
        diff.missing_in_baseline = missing_in_baseline;
    end

    % Compare common fields
    common_fields = intersect(data_fields, baseline_fields);
    diff.field_diffs = {};

    for i = 1:length(common_fields)
        field = common_fields{i};
        data_val = data.(field);
        baseline_val = baseline.(field);

        if ~compare_values(data_val, baseline_val, tolerance)
            match = false;
            diff.field_diffs{end+1} = field;
        end
    end
end

function match = compare_values(val1, val2, tolerance)
    % Compare two values with tolerance for numeric types
    match = true;

    if ~strcmp(class(val1), class(val2))
        match = false;
        return;
    end

    if isnumeric(val1) && isnumeric(val2)
        if ~isequal(size(val1), size(val2))
            match = false;
            return;
        end
        if any(abs(val1(:) - val2(:)) > tolerance, 'all')
            match = false;
        end
    elseif ischar(val1) && ischar(val2)
        match = strcmp(val1, val2);
    elseif isstruct(val1) && isstruct(val2)
        % Recursive comparison for structs
        fields1 = fieldnames(val1);
        fields2 = fieldnames(val2);
        if ~isequal(sort(fields1), sort(fields2))
            match = false;
            return;
        end
        for i = 1:length(fields1)
            if ~compare_values(val1.(fields1{i}), val2.(fields1{i}), tolerance)
                match = false;
                return;
            end
        end
    elseif iscell(val1) && iscell(val2)
        if ~isequal(size(val1), size(val2))
            match = false;
            return;
        end
        for i = 1:numel(val1)
            if ~compare_values(val1{i}, val2{i}, tolerance)
                match = false;
                return;
            end
        end
    else
        match = isequal(val1, val2);
    end
end

function results = list_baselines(baselines_dir)
    % List all baseline files
    results = struct();
    results.files = {};

    if ~exist(baselines_dir, 'dir')
        fprintf('No baselines directory found: %s\n', baselines_dir);
        return;
    end

    files = dir(fullfile(baselines_dir, '**', '*.mat'));
    fprintf('Baseline files (%d total):\n', length(files));

    for i = 1:length(files)
        rel_path = strrep(files(i).folder, baselines_dir, '');
        if ~isempty(rel_path) && (rel_path(1) == filesep)
            rel_path = rel_path(2:end);
        end
        if isempty(rel_path)
            full_name = files(i).name;
        else
            full_name = fullfile(rel_path, files(i).name);
        end
        results.files{end+1} = full_name;
        fprintf('  %s\n', full_name);
    end
end

function results = check_status(data_dir, baselines_dir)
    % Check status of baselines vs current data
    results = struct();
    results.has_baselines = exist(baselines_dir, 'dir') && ...
        ~isempty(dir(fullfile(baselines_dir, '**', '*.mat')));
    results.has_data = exist(data_dir, 'dir') && ...
        ~isempty(dir(fullfile(data_dir, '**', '*.mat')));

    fprintf('Baseline Status:\n');
    fprintf('  Baselines directory exists: %d\n', exist(baselines_dir, 'dir') > 0);
    fprintf('  Data directory exists: %d\n', exist(data_dir, 'dir') > 0);

    if results.has_baselines
        baseline_files = dir(fullfile(baselines_dir, '**', '*.mat'));
        fprintf('  Baseline files: %d\n', length(baseline_files));
    else
        fprintf('  Baseline files: 0 (run manage_baselines(''save'') to create)\n');
    end

    if results.has_data
        data_files = dir(fullfile(data_dir, '**', '*.mat'));
        fprintf('  Data files: %d\n', length(data_files));
    else
        fprintf('  Data files: 0 (run tests first)\n');
    end
end

function results = clean_data(data_dir)
    % Clean current test data (not baselines)
    results = struct();
    results.deleted = {};

    if ~exist(data_dir, 'dir')
        fprintf('Data directory does not exist: %s\n', data_dir);
        return;
    end

    files = dir(fullfile(data_dir, '**', '*.mat'));
    for i = 1:length(files)
        file_path = fullfile(files(i).folder, files(i).name);
        delete(file_path);
        results.deleted{end+1} = files(i).name;
    end

    fprintf('Deleted %d data files.\n', length(results.deleted));
end
