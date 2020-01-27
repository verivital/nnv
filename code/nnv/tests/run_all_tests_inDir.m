% execute all tests in a given directory, pausing for t_pause between each (as many show plots,
% etc.), assuming test files start or end with "test_" or "_test"
function [disabledTests, i_d] = run_all_tests_inDir(dirname, outputDirname)
    global disabledTests i_d;

    dirname

    cd(char(dirname)); % char necessary for linux/codeocean
    fprintf('\n\nEXECUTING TESTS IN: %s\n', pwd);
    t_pause = 0.01;
    
    % don't run these tests (throw exceptions)
    % TODO: fix these exceptions and tests, we cannot have tests 
    % with exceptions for codeocean
    offTests = {'test_box_getVertices_c_gen_mexexa64', 'test_box_getVertices_c_gen_rtwk', 'test_FFNNS_reach_star', 'test_NonLinearODE_evaluate', 'test_NonLinearODE_reach_zono', ...
            'test_star_plotBoxes_3D', 'test_PosLin_stepReach', 'test_ReLU_reach_approx', 'test_NonLinearODE_stepReachStar', 'test_LayerS_reach_satlin',  'test_LayerS_reach_poslin', 'test_CNN_parse'};

    % to manually get all files in reasonable format for setting up tests
    % fprintf("%s\n",strtrim(strrep(strrep(string(ls),'"',''),'.m','')))

    % windows
    if ispc
    	test_files = strtrim(strrep(strrep(string(ls),'"',''),'.m',''));
    else % unix/codeocean
    	test_files = strtrim(strsplit(strrep(strrep(string(ls),'"',''),'.m','')));
    end
    
    test_files

    % iterate over all file names starting with "test_" or ending with "_test" and execute them
    for i_f = 1:length(test_files)
        close all ; clearAllExceptVars({t_pause,test_files,i_f,disabledTests,i_d}); % save a few variables between tests, otherwise clear
        f = test_files(i_f);
        
        if (strncmpi(f, 'test_',5) == 1 || strncmpi(reverse(f), reverse('_test'),5) == 1) && ~contains(f,'.') && sum(contains(offTests, f)) == 0
            fprintf('running test %s\n', f);
            try
                run(f);
            catch e
                fprintf('\n\nERROR running test: %s\n', f);
                disabledTests{i_d} = f; % add to list of test names that should be disabled (for codeocean)
                %fprintf('%s\n\n', e); % show error message
                e
                i_d = i_d + 1;
            end
        end
        pause(t_pause);
        
        % save all figures to files
        figHandles = findall(groot, 'Type', 'figure');

        for i_fh = 1 : length(figHandles)
            fh = figHandles(i_fh);
            filename = strcat(outputDirname, 'results_fig_', f, num2str(i_fh), '.png');
            saveas(fh, filename);
        end
        close all;
    end
    
    close all ; clearAllExceptVars({t_pause,disabledTests,i_d});
    cd('..');
end