% execute all tests recursively on subdirectories from this directory
%
% example run:
% clear all ; disabledTests = {}; i_d = 1; global disabledTests i_d; close all ; clc ; run_all_tests([pwd, filesep, '..', filesep, '..', filesep, '..', filesep, 'results', filesep], disabledTests, i_d)
function all_results=run_all_tests(outputDirname, disabledTests, i_d)
    global disabledTests i_d;

    % windows
    if ispc
    	dirNames = strtrim(strrep(string(ls),'"',''));
    else % linux/codeocean
    	dirNames = strtrim(strsplit(strrep(string(ls),'"','')));
    end
    
    ls
    dirNames
    all_results=[]
    for i_d = 1 : length(dirNames)
        d = dirNames(i_d);
        if isdir(d) && ~(strcmp(d,'.') || strcmp(d,'..'))
	  d
            run_all_tests_inDir(d, outputDirname);
            
            d
            cd(char(d)); % char necessary for linux/codeocean
            [results, g1, g2]=run_all_tests(outputDirname, disabledTests, i_d); % recursive call
	    all_results=[all_results, results];
            cd('..');
        end
    end
end
