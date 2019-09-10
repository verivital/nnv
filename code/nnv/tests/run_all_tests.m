% execute all tests recursively on subdirectories from this directory
%
% example run:
% clear all ; disabledTests = {}; i_d = 1; global disabledTests i_d; close all ; clc ; run_all_tests([pwd, filesep, '..', filesep, '..', filesep, '..', filesep, 'results', filesep], disabledTests, i_d)
function run_all_tests(outputDirname, disabledTests, i_d)
    global disabledTests i_d;

    % windows
    if ispc
    	dirNames = strtrim(strrep(string(ls),'"',''));
    else % linux/codeocean
    	dirNames = strtrim(strsplit(strrep(string(ls),'"','')));
    end
    
    ls
    dirNames
    
    for i_d = 1 : length(dirNames)
        d = dirNames(i_d);
        if isdir(d) && ~(strcmp(d,'.') || strcmp(d,'..'))
            run_all_tests_inDir(d, outputDirname);
            
            d
            cd(char(d)); % char necessary for linux/codeocean
            run_all_tests(outputDirname, disabledTests, i_d); % recursive call
            cd('..');
        end
    end
end