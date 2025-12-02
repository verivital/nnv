% run_all_tests.m
% Master test runner for all NNV tests
% Unified entry point for soundness, regression, and utility tests
%
% Usage:
%   results = run_all_tests           % Run all tests with concise output
%   results = run_all_tests('verbose')  % Run with detailed output
%   results = run_all_tests('quick')    % Run only soundness/regression/utils
%
% Returns:
%   results - MATLAB test results array
%   Exit code 0 if all tests pass, 1 if any fail (for CI/CD)

function results = run_all_tests(varargin)
    verbose = nargin > 0 && strcmpi(varargin{1}, 'verbose');
    quick = nargin > 0 && strcmpi(varargin{1}, 'quick');

    start_time = tic;

    fprintf('\n');
    fprintf('=======================================================\n');
    fprintf('           NNV COMPLETE TEST SUITE                     \n');
    fprintf('=======================================================\n');
    fprintf('Started: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

    % Get the tests directory
    test_dir = fileparts(mfilename('fullpath'));

    % Define test categories (ordered by priority)
    if quick
        % Quick mode: only core test suites
        categories = {
            'soundness', 'Soundness Tests'
            'regression', 'Regression Tests'
            'utils',     'Utility Function Tests'
        };
        fprintf('Running in QUICK mode (soundness + regression + utils only)\n\n');
    else
        % Full mode: all test categories
        categories = {
            'soundness', 'Soundness Tests'
            'regression', 'Regression Tests'
            'utils',     'Utility Function Tests'
            'nn',        'Neural Network Tests'
            'nncs',      'Neural Network Control Systems Tests'
            'set',       'Set Representation Tests'
            'glpk',      'GLPK Tests'
            'io',        'I/O Tests'
            'tutorial',  'Tutorial Tests'
        };
    end

    all_results = [];
    category_summary = {};

    for i = 1:size(categories, 1)
        folder = categories{i, 1};
        label = categories{i, 2};
        folder_path = fullfile(test_dir, folder);

        if ~isfolder(folder_path)
            fprintf('Skipping %s (folder not found)\n', label);
            continue;
        end

        fprintf('-------------------------------------------------------\n');
        fprintf('Running: %s\n', label);
        fprintf('Path: %s\n', folder_path);
        fprintf('-------------------------------------------------------\n');

        try
            % Run tests in this category
            if verbose
                cat_results = runtests(folder_path, 'IncludeSubfolders', true);
            else
                cat_results = runtests(folder_path, 'IncludeSubfolders', true, 'OutputDetail', 'Concise');
            end

            if ~isempty(cat_results)
                all_results = [all_results, cat_results];

                passed = sum([cat_results.Passed]);
                failed = sum([cat_results.Failed]);
                total = length(cat_results);

                category_summary{end+1, 1} = label;
                category_summary{end, 2} = total;
                category_summary{end, 3} = passed;
                category_summary{end, 4} = failed;

                fprintf('%s: %d/%d passed\n\n', label, passed, total);
            else
                fprintf('%s: No tests found\n\n', label);
            end
        catch ME
            fprintf('ERROR running %s: %s\n\n', label, ME.message);
            category_summary{end+1, 1} = label;
            category_summary{end, 2} = 0;
            category_summary{end, 3} = 0;
            category_summary{end, 4} = 0;
        end
    end

    % Print overall summary
    fprintf('\n');
    fprintf('=======================================================\n');
    fprintf('                 OVERALL SUMMARY                       \n');
    fprintf('=======================================================\n\n');

    if ~isempty(all_results)
        total = length(all_results);
        passed = sum([all_results.Passed]);
        failed = sum([all_results.Failed]);
        incomplete = sum([all_results.Incomplete]);

        fprintf('%-40s %8s %8s %8s\n', 'Category', 'Total', 'Passed', 'Failed');
        fprintf('%-40s %8s %8s %8s\n', repmat('-', 1, 40), '------', '------', '------');

        for i = 1:size(category_summary, 1)
            fprintf('%-40s %8d %8d %8d\n', ...
                category_summary{i, 1}, ...
                category_summary{i, 2}, ...
                category_summary{i, 3}, ...
                category_summary{i, 4});
        end

        fprintf('%-40s %8s %8s %8s\n', repmat('-', 1, 40), '------', '------', '------');
        fprintf('%-40s %8d %8d %8d\n', 'TOTAL', total, passed, failed);
        fprintf('\n');

        if failed == 0 && incomplete == 0
            fprintf('*** ALL %d TESTS PASSED! ***\n', total);
        else
            fprintf('*** %d TESTS FAILED, %d INCOMPLETE ***\n', failed, incomplete);

            % List failed tests
            fprintf('\nFailed tests:\n');
            for i = 1:length(all_results)
                if all_results(i).Failed
                    fprintf('  - %s\n', all_results(i).Name);
                end
            end
        end
    else
        fprintf('No tests were run.\n');
    end

    % Print elapsed time
    elapsed = toc(start_time);
    fprintf('\nCompleted: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf('Elapsed time: %.1f seconds (%.1f minutes)\n', elapsed, elapsed/60);
    fprintf('\n');

    results = all_results;
end
