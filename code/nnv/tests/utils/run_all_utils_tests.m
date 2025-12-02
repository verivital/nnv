% run_all_utils_tests.m
% Runner script for NNV utility function tests
% Usage: run('run_all_utils_tests.m')

fprintf('=== Running NNV Utility Function Tests ===\n\n');

% Get the directory containing this script
test_dir = fileparts(mfilename('fullpath'));
cd(test_dir);

% Find all test files
test_files = dir('test_*.m');

total_tests = 0;
passed_tests = 0;
failed_tests = 0;
failed_names = {};

for i = 1:length(test_files)
    test_name = test_files(i).name(1:end-2);  % Remove .m extension
    fprintf('Running %s...\n', test_name);

    try
        results = runtests(test_name);

        % Count results
        for j = 1:length(results)
            total_tests = total_tests + 1;
            if results(j).Passed
                passed_tests = passed_tests + 1;
            else
                failed_tests = failed_tests + 1;
                failed_names{end+1} = [test_name, '/', results(j).Name];
            end
        end
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        failed_tests = failed_tests + 1;
        failed_names{end+1} = test_name;
    end
end

fprintf('\n=== Summary ===\n');
fprintf('Total: %d, Passed: %d, Failed: %d\n', total_tests, passed_tests, failed_tests);

if failed_tests > 0
    fprintf('\nFailed tests:\n');
    for i = 1:length(failed_names)
        fprintf('  - %s\n', failed_names{i});
    end
else
    fprintf('All utility tests passed!\n');
end
