function out = is_github_actions()
    nnvpath = nnvroot();
    file_path = [nnvpath, 'github_actions.txt'];
    % for debugging
    disp(nnvpath);
    disp(file_path);
    allfiles = dir(nnvpath);
    tableFiles = struct2table(allfiles);
    disp(tableFiles);
    if isfile(file_path)
        out = 1;
        % 'github actions detected'
    else
        out = 0;
        % 'github actions not detected'
    end
end

