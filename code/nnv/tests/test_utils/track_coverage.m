function results = track_coverage(varargin)
    % TRACK_COVERAGE - Run NNV tests with code coverage tracking
    %
    % This function runs NNV tests while collecting code coverage data
    % using MATLAB's built-in profiler and coverage tools.
    %
    % Usage:
    %   results = track_coverage()                     % Run all tests with coverage
    %   results = track_coverage('quick')              % Quick mode (fewer tests)
    %   results = track_coverage('verbose')            % Verbose output
    %   results = track_coverage('report', true)       % Generate HTML report
    %   results = track_coverage('folder', 'tests/nn') % Specific test folder
    %
    % Options:
    %   'mode'    - 'full', 'quick', or 'verbose' (default: 'full')
    %   'report'  - Generate HTML coverage report (default: false)
    %   'folder'  - Test folder to run (default: 'tests')
    %   'source'  - Source folders to track (default: {'engine'})
    %
    % Returns:
    %   results - Struct with fields:
    %       .test_results    - MATLAB test results
    %       .coverage        - Coverage statistics struct
    %       .summary         - Text summary string
    %
    % Example:
    %   % Run quick tests with coverage report
    %   r = track_coverage('quick', 'report', true);
    %   fprintf('Coverage: %.1f%%\n', r.coverage.percent);
    %
    % See also: run_all_tests, profile, runtests

    % Parse arguments
    p = inputParser;
    addOptional(p, 'mode', 'full', @(x) ischar(x) && ismember(x, {'full', 'quick', 'verbose'}));
    addParameter(p, 'report', false, @islogical);
    addParameter(p, 'folder', '', @ischar);
    addParameter(p, 'source', {'engine'}, @iscell);
    parse(p, varargin{:});
    opts = p.Results;

    % Get NNV root directory
    nnv_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

    fprintf('\n');
    fprintf('=======================================================\n');
    fprintf('         NNV TEST COVERAGE ANALYSIS                    \n');
    fprintf('=======================================================\n');
    fprintf('Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('NNV Root: %s\n', nnv_root);
    fprintf('Mode: %s\n\n', opts.mode);

    %% Identify source files to track
    fprintf('Identifying source files to track...\n');
    source_files = {};

    for i = 1:length(opts.source)
        source_dir = fullfile(nnv_root, opts.source{i});
        if isfolder(source_dir)
            % Get all .m files in this source directory
            m_files = dir(fullfile(source_dir, '**', '*.m'));
            for j = 1:length(m_files)
                full_path = fullfile(m_files(j).folder, m_files(j).name);
                % Exclude test files and external code
                if ~contains(full_path, 'test') && ...
                   ~contains(full_path, 'cora') && ...
                   ~contains(full_path, 'tbxmanager')
                    source_files{end+1} = full_path;
                end
            end
        end
    end

    fprintf('  Found %d source files to track\n\n', length(source_files));

    %% Set up coverage collection using profiler
    fprintf('Setting up coverage collection...\n');

    % Clear any existing profile data
    profile off;
    profile clear;

    % Start profiler with coverage tracking
    profile on -history;

    %% Run tests
    fprintf('Running tests with coverage enabled...\n\n');
    test_start = tic;

    try
        test_dir = fullfile(nnv_root, 'tests');

        if ~isempty(opts.folder)
            test_dir = fullfile(nnv_root, opts.folder);
        end

        % Configure test runner based on mode
        switch opts.mode
            case 'quick'
                % Run only soundness and regression tests
                soundness_dir = fullfile(nnv_root, 'tests', 'soundness');
                regression_dir = fullfile(nnv_root, 'tests', 'regression');

                test_results1 = runtests(soundness_dir, 'IncludeSubfolders', true, 'OutputDetail', 'Concise');
                test_results2 = runtests(regression_dir, 'IncludeSubfolders', true, 'OutputDetail', 'Concise');
                test_results = [test_results1, test_results2];

            case 'verbose'
                test_results = runtests(test_dir, 'IncludeSubfolders', true);

            otherwise  % 'full'
                test_results = runtests(test_dir, 'IncludeSubfolders', true, 'OutputDetail', 'Concise');
        end
    catch ME
        fprintf('Error running tests: %s\n', ME.message);
        test_results = [];
    end

    test_time = toc(test_start);

    %% Stop profiler and collect data
    profile off;
    prof_data = profile('info');

    %% Analyze coverage
    fprintf('\n');
    fprintf('=======================================================\n');
    fprintf('         COVERAGE ANALYSIS                             \n');
    fprintf('=======================================================\n\n');

    coverage = analyze_coverage(prof_data, source_files, nnv_root);

    %% Print coverage summary
    fprintf('Coverage by Directory:\n');
    fprintf('%-50s %8s %8s %10s\n', 'Directory', 'Files', 'Executed', 'Coverage');
    fprintf('%-50s %8s %8s %10s\n', repmat('-', 1, 50), '------', '--------', '--------');

    dir_names = fieldnames(coverage.by_directory);
    for i = 1:length(dir_names)
        dir_info = coverage.by_directory.(dir_names{i});
        fprintf('%-50s %8d %8d %9.1f%%\n', ...
            dir_names{i}, dir_info.total, dir_info.executed, dir_info.percent);
    end

    fprintf('%-50s %8s %8s %10s\n', repmat('-', 1, 50), '------', '--------', '--------');
    fprintf('%-50s %8d %8d %9.1f%%\n', 'TOTAL', ...
        coverage.total_files, coverage.executed_files, coverage.percent);
    fprintf('\n');

    %% List uncovered files (most important ones)
    if ~isempty(coverage.uncovered)
        fprintf('Key Uncovered Files (first 20):\n');
        show_count = min(20, length(coverage.uncovered));
        for i = 1:show_count
            % Make path relative
            rel_path = strrep(coverage.uncovered{i}, [nnv_root filesep], '');
            fprintf('  - %s\n', rel_path);
        end
        if length(coverage.uncovered) > 20
            fprintf('  ... and %d more\n', length(coverage.uncovered) - 20);
        end
        fprintf('\n');
    end

    %% Generate HTML report if requested
    if opts.report
        report_dir = fullfile(nnv_root, 'tests', 'coverage_report');
        fprintf('Generating HTML coverage report...\n');

        try
            if ~isfolder(report_dir)
                mkdir(report_dir);
            end

            % Use MATLAB's profiler report
            profreport(report_dir);
            fprintf('  Report saved to: %s\n\n', report_dir);
        catch ME
            fprintf('  Could not generate HTML report: %s\n\n', ME.message);
        end
    end

    %% Build results struct
    results = struct();
    results.test_results = test_results;
    results.coverage = coverage;
    results.elapsed_time = test_time;

    % Build summary string
    if ~isempty(test_results)
        n_passed = sum([test_results.Passed]);
        n_failed = sum([test_results.Failed]);
        n_total = length(test_results);
    else
        n_passed = 0;
        n_failed = 0;
        n_total = 0;
    end

    results.summary = sprintf('Tests: %d/%d passed | Coverage: %.1f%% (%d/%d files)', ...
        n_passed, n_total, coverage.percent, coverage.executed_files, coverage.total_files);

    %% Final summary
    fprintf('=======================================================\n');
    fprintf('                 FINAL SUMMARY                         \n');
    fprintf('=======================================================\n');

    if ~isempty(test_results)
        fprintf('Tests:    %d passed, %d failed (%.1f%% pass rate)\n', ...
            n_passed, n_failed, 100 * n_passed / max(1, n_total));
    end

    fprintf('Coverage: %.1f%% (%d of %d source files executed)\n', ...
        coverage.percent, coverage.executed_files, coverage.total_files);
    fprintf('Time:     %.1f seconds (%.1f minutes)\n', test_time, test_time/60);
    fprintf('=======================================================\n\n');
end


function coverage = analyze_coverage(prof_data, source_files, nnv_root)
    % Analyze profiler data to determine code coverage

    % Get list of all files that were executed
    executed_files = {};
    if ~isempty(prof_data.FunctionTable)
        executed_files = unique({prof_data.FunctionTable.FileName});
        % Filter to only .m files
        executed_files = executed_files(cellfun(@(x) endsWith(x, '.m'), executed_files));
    end

    % Normalize paths for comparison
    source_files_norm = cellfun(@(x) lower(x), source_files, 'UniformOutput', false);
    executed_files_norm = cellfun(@(x) lower(x), executed_files, 'UniformOutput', false);

    % Count executed source files
    executed_count = 0;
    uncovered = {};

    for i = 1:length(source_files)
        if any(strcmp(source_files_norm{i}, executed_files_norm))
            executed_count = executed_count + 1;
        else
            uncovered{end+1} = source_files{i};
        end
    end

    % Calculate overall coverage
    coverage = struct();
    coverage.total_files = length(source_files);
    coverage.executed_files = executed_count;
    coverage.percent = 100 * executed_count / max(1, length(source_files));
    coverage.uncovered = uncovered;

    % Calculate coverage by directory
    coverage.by_directory = struct();

    % Define key directories to track
    key_dirs = {'nn', 'nncs', 'set', 'utils', 'hyst', 'cnn'};

    for i = 1:length(key_dirs)
        dir_name = key_dirs{i};
        dir_pattern = [filesep dir_name filesep];

        % Count files in this directory
        in_dir = cellfun(@(x) contains(lower(x), lower(dir_pattern)), source_files);
        executed_in_dir = cellfun(@(x) contains(lower(x), lower(dir_pattern)), source_files) & ...
                          cellfun(@(x) any(strcmp(lower(x), executed_files_norm)), source_files);

        dir_total = sum(in_dir);
        dir_executed = sum(executed_in_dir);

        if dir_total > 0
            % Use valid field name
            field_name = matlab.lang.makeValidName(dir_name);
            coverage.by_directory.(field_name) = struct(...
                'total', dir_total, ...
                'executed', dir_executed, ...
                'percent', 100 * dir_executed / dir_total);
        end
    end
end
