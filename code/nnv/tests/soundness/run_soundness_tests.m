function results = run_soundness_tests(varargin)
    % RUN_SOUNDNESS_TESTS Run comprehensive soundness tests for SLM layers
    %
    % Usage:
    %   results = run_soundness_tests()           % Run all soundness tests
    %   results = run_soundness_tests('quick')    % Run subset of tests
    %   results = run_soundness_tests('all')      % Run all soundness tests
    %
    % Author: NNV Team
    % Date: December 2025

    % Parse mode
    if nargin > 0
        mode = lower(varargin{1});
    else
        mode = 'all';
    end

    % Define test files
    core_tests = {
        'test_SiLULayer_soundness'
        'test_SwiGLULayer_soundness'
        'test_RMSNormLayer_soundness'
        'test_AdditionLayer_soundness'
    };

    % Get the directory containing soundness tests
    test_dir = fileparts(mfilename('fullpath'));

    % Store current directory and change to test directory
    orig_dir = pwd;
    cd(test_dir);

    try
        fprintf('\n');
        fprintf('==========================================\n');
        fprintf('     NNV Soundness Test Suite            \n');
        fprintf('==========================================\n');
        fprintf('Mode: %s\n', mode);
        fprintf('Time: %s\n', datestr(now));
        fprintf('==========================================\n\n');

        switch mode
            case 'quick'
                fprintf('[Quick Mode] Running core soundness tests...\n\n');
                tests_to_run = core_tests(1:2);  % Just SiLU and SwiGLU

            otherwise  % 'all' or default
                fprintf('[All Mode] Running all soundness tests...\n\n');
                tests_to_run = core_tests;
        end

        % Check which tests exist
        existing_tests = {};
        missing_tests = {};
        for i = 1:length(tests_to_run)
            test_file = [tests_to_run{i}, '.m'];
            if exist(test_file, 'file')
                existing_tests{end+1} = tests_to_run{i};
            else
                missing_tests{end+1} = tests_to_run{i};
            end
        end

        % Report missing tests
        if ~isempty(missing_tests)
            fprintf('Warning: %d test files not found:\n', length(missing_tests));
            for i = 1:length(missing_tests)
                fprintf('  - %s.m\n', missing_tests{i});
            end
            fprintf('\n');
        end

        % Run tests
        if isempty(existing_tests)
            fprintf('No test files found!\n');
            results = [];
        else
            fprintf('Running %d soundness test files:\n', length(existing_tests));
            for i = 1:length(existing_tests)
                fprintf('  %d. %s\n', i, existing_tests{i});
            end
            fprintf('\n');

            % Run the tests
            results = runtests(existing_tests);

            % Summary
            fprintf('\n');
            fprintf('==========================================\n');
            fprintf('         SOUNDNESS TEST SUMMARY          \n');
            fprintf('==========================================\n');

            n_passed = sum([results.Passed]);
            n_failed = sum([results.Failed]);
            n_incomplete = sum([results.Incomplete]);
            n_total = length(results);

            fprintf('Total Tests: %d\n', n_total);
            fprintf('Passed:      %d (%.1f%%)\n', n_passed, 100*n_passed/n_total);
            fprintf('Failed:      %d\n', n_failed);
            fprintf('Incomplete:  %d\n', n_incomplete);
            fprintf('==========================================\n');

            if n_failed == 0 && n_incomplete == 0
                fprintf('\n*** ALL SOUNDNESS TESTS PASSED ***\n');
                fprintf('All reachability analyses are SOUND.\n\n');
            else
                fprintf('\n*** SOME SOUNDNESS TESTS FAILED ***\n');
                fprintf('Review failures - potential soundness violations!\n\n');
            end
        end

    catch ME
        cd(orig_dir);
        rethrow(ME);
    end

    % Restore original directory
    cd(orig_dir);
end
