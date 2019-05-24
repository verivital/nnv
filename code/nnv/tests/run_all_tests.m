% execute all tests recursively on subdirectories from this directory
function run_all_tests(outputDirname, disabledTests, i_d)
    global disabledTests i_d;
    ' starting tests'
    fprintf('starting all tests in root\n');
    dirNames = strtrim(strsplit(strrep(string(ls),'"','')));
    
    strtrim(strsplit(strrep(string(ls),'"','')))
    ls
    dirNames
    for i_d = 1 : length(dirNames)
        d = dirNames(i_d);
        d
        dirNames(i_d)
        i_d
        if isdir(d) && ~(strcmp(d,'.') || strcmp(d,'..'))
            run_all_tests_inDir(d, outputDirname);
            
            cd(char(d));
            run_all_tests(outputDirname, disabledTests, i_d); % recursive call
            cd('..');
        end
    end
end