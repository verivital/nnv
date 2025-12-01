% run_all_soundness_tests.m
% Runs all soundness tests in this folder and reports summary
% Usage: run_all_soundness_tests

function results = run_all_soundness_tests()
    fprintf('\n=== NNV SOUNDNESS TEST SUITE ===\n\n');

    % Get the directory containing this script
    test_dir = fileparts(mfilename('fullpath'));

    % Find all test files (test_soundness_*.m)
    test_files = dir(fullfile(test_dir, 'test_soundness_*.m'));

    if isempty(test_files)
        fprintf('No soundness test files found in: %s\n', test_dir);
        results = [];
        return;
    end

    fprintf('Found %d soundness test files:\n', length(test_files));
    for i = 1:length(test_files)
        fprintf('  %d. %s\n', i, test_files(i).name);
    end
    fprintf('\n');

    % Run all tests
    fprintf('Running tests...\n\n');
    results = runtests(test_dir, 'IncludeSubfolders', false);

    % Print summary
    fprintf('\n=== SOUNDNESS TEST SUMMARY ===\n');

    total = length(results);
    passed = sum([results.Passed]);
    failed = sum([results.Failed]);
    incomplete = sum([results.Incomplete]);

    fprintf('Total:      %d\n', total);
    fprintf('Passed:     %d\n', passed);
    fprintf('Failed:     %d\n', failed);
    fprintf('Incomplete: %d\n', incomplete);
    fprintf('\n');

    if failed == 0 && incomplete == 0
        fprintf('ALL SOUNDNESS TESTS PASSED!\n');
    else
        fprintf('SOME TESTS FAILED OR INCOMPLETE.\n');
        fprintf('\nFailed/Incomplete tests:\n');
        for i = 1:length(results)
            if results(i).Failed || results(i).Incomplete
                fprintf('  - %s\n', results(i).Name);
            end
        end
    end
    fprintf('\n');
end
