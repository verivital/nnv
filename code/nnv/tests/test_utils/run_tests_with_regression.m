function [test_results, regression_results] = run_tests_with_regression(varargin)
    % RUN_TESTS_WITH_REGRESSION - Run NNV tests with regression detection
    %
    % This script runs the NNV test suite and compares outputs against
    % saved baselines to detect regressions.
    %
    % Usage:
    %   run_tests_with_regression()                    - Run all tests
    %   run_tests_with_regression('tests/set')         - Run specific test folder
    %   run_tests_with_regression('compare', true)     - Force baseline comparison
    %   run_tests_with_regression('compare', false)    - Skip baseline comparison
    %   run_tests_with_regression('verbose', true)     - Verbose output
    %
    % Options:
    %   'folder'    - Test folder to run (default: 'tests')
    %   'compare'   - Compare against baselines (default: from config)
    %   'verbose'   - Verbose output (default: false)
    %   'parallel'  - Use parallel execution (default: false)
    %
    % Returns:
    %   test_results      - MATLAB test results object
    %   regression_results - Baseline comparison results struct
    %
    % Example (CI/CD):
    %   setenv('NNV_TEST_COMPARE_BASELINES', '1');
    %   [results, regressions] = run_tests_with_regression();
    %   if ~isempty(regressions.regressions)
    %       exit(1);  % Fail CI
    %   end

    % Parse arguments
    p = inputParser;
    addOptional(p, 'folder', 'tests', @ischar);
    addParameter(p, 'compare', [], @(x) islogical(x) || isempty(x));
    addParameter(p, 'verbose', false, @islogical);
    addParameter(p, 'parallel', false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;

    % Get test configuration
    config = get_test_config();

    % Determine if we should compare baselines
    if isempty(opts.compare)
        do_compare = config.compare_baselines;
    else
        do_compare = opts.compare;
    end

    fprintf('========================================\n');
    fprintf('NNV Test Suite with Regression Detection\n');
    fprintf('========================================\n');
    fprintf('Test folder: %s\n', opts.folder);
    fprintf('Baseline comparison: %s\n', mat2str(do_compare));
    fprintf('Verbose: %s\n', mat2str(opts.verbose));
    fprintf('Timestamp: %s\n', datestr(now));
    fprintf('========================================\n\n');

    % Run the tests
    fprintf('Running tests...\n');
    test_start = tic;

    try
        if opts.parallel
            test_results = runtests(opts.folder, 'IncludeSubfolders', true, ...
                'UseParallel', true);
        else
            test_results = runtests(opts.folder, 'IncludeSubfolders', true);
        end
    catch ME
        fprintf('Error running tests: %s\n', ME.message);
        test_results = [];
        regression_results = struct('error', ME.message);
        return;
    end

    test_time = toc(test_start);
    fprintf('\nTests completed in %.1f seconds.\n', test_time);

    % Summarize test results
    if ~isempty(test_results)
        n_passed = sum([test_results.Passed]);
        n_failed = sum([test_results.Failed]);
        n_incomplete = sum([test_results.Incomplete]);

        fprintf('\nTest Results Summary:\n');
        fprintf('  Passed: %d\n', n_passed);
        fprintf('  Failed: %d\n', n_failed);
        fprintf('  Incomplete: %d\n', n_incomplete);
        fprintf('  Total: %d\n', length(test_results));
    end

    % Compare against baselines if enabled
    regression_results = struct();
    regression_results.compared = do_compare;
    regression_results.regressions = {};

    if do_compare
        fprintf('\n========================================\n');
        fprintf('Comparing against baselines...\n');
        fprintf('========================================\n');

        regression_results = manage_baselines('compare', 'verbose', opts.verbose);

        if ~isempty(regression_results.regressions)
            fprintf('\n*** REGRESSION DETECTED ***\n');
            fprintf('The following tests produced different results than baseline:\n');
            for i = 1:length(regression_results.regressions)
                fprintf('  - %s\n', regression_results.regressions{i}.file);
            end

            if config.fail_on_regression
                fprintf('\nFailing due to regression (config.fail_on_regression = true)\n');
                fprintf('To update baselines, run: manage_baselines(''save'')\n');
            end
        else
            fprintf('\nNo regressions detected. All outputs match baselines.\n');
        end
    else
        fprintf('\nBaseline comparison skipped.\n');
        fprintf('To enable: set compare_baselines=true in get_test_config.m\n');
        fprintf('Or run: run_tests_with_regression(''compare'', true)\n');
    end

    % Final summary
    fprintf('\n========================================\n');
    fprintf('Final Summary\n');
    fprintf('========================================\n');

    all_passed = true;

    if ~isempty(test_results)
        if sum([test_results.Failed]) > 0
            fprintf('TEST FAILURES: %d tests failed\n', sum([test_results.Failed]));
            all_passed = false;
        else
            fprintf('All %d tests passed.\n', sum([test_results.Passed]));
        end
    end

    if do_compare && ~isempty(regression_results.regressions)
        fprintf('REGRESSIONS: %d baseline mismatches\n', length(regression_results.regressions));
        all_passed = false;
    end

    if all_passed
        fprintf('\n*** ALL CHECKS PASSED ***\n');
    else
        fprintf('\n*** CHECKS FAILED ***\n');
    end

    fprintf('========================================\n');
end
