%% run_all_regression_tests
% Master script to run all regression tests
% Usage: run('run_all_regression_tests.m')

fprintf('=== NNV Regression Tests ===\n');
fprintf('Started at: %s\n\n', datestr(now));

% Store current directory
start_dir = pwd;

% Change to regression test directory
test_dir = fileparts(mfilename('fullpath'));
cd(test_dir);

% Run all regression tests
results = runtests(pwd, 'IncludeSubfolders', false);

% Return to original directory
cd(start_dir);

% Display summary
fprintf('\n=== Regression Test Summary ===\n');

passed = sum([results.Passed]);
failed = sum([results.Failed]);
incomplete = sum([results.Incomplete]);
total = length(results);

fprintf('Total:      %d\n', total);
fprintf('Passed:     %d\n', passed);
fprintf('Failed:     %d\n', failed);
fprintf('Incomplete: %d\n', incomplete);

if failed > 0
    fprintf('\n=== Failed Tests ===\n');
    for i = 1:length(results)
        if results(i).Failed
            fprintf('  - %s\n', results(i).Name);
        end
    end
end

fprintf('\nCompleted at: %s\n', datestr(now));

% Return results for programmatic use
if nargout > 0
    varargout{1} = results;
end
