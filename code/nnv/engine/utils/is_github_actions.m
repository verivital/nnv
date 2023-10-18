function out = is_github_actions()
    if isfolder('/github_actions-true')
        out = 1;
        % 'github actions detected'
    else
        out = 0;
        % 'github actions not detected'
    end
end

