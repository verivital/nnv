% execute all tests in a given directory, pausing for t_pause between each (as many show plots,
% etc.), assuming test files start or end with "test_" or "_test"
function run_all_tests_inDir(dirname)
    cd(dirname);
    fprintf('\n\nEXECUTING TESTS IN: %s\n', pwd);
    t_pause = 0.1;

    % to manually get all files in reasonable format for setting up tests
    % fprintf("%s\n",strtrim(strrep(strrep(string(ls),'"',''),'.m','')))

    test_files = strtrim(strrep(strrep(string(ls),'"',''),'.m',''));

    % iterate over all file names starting with "test_" or ending with "_test" and execute them
    for i_f = 1:length(test_files)
        close all ; clearAllExceptVars({t_pause,test_files,i_f}); % save a few variables between tests, otherwise clear
        f = test_files(i_f);
        if (strncmpi(f, 'test_',5) == 1 || strncmpi(reverse(f), reverse('_test'),5) == 1) && ~contains(f,'.')
            fprintf('running test %s\n', f);
            try
                run(f);
            catch e
                fprintf('\n\nERROR running test: %s\n', f);
                e % show error message
            end
        end
        pause(t_pause);
    end

    close all ; clear all;
    cd('..');
end