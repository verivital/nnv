function randomEx_table()

    % Generate time tables for paper submission
    files = dir('results/*.mat'); % order alphabetically [L,M,S,XL,XS,XXL]
    res = [];
    
    % Load files
    for i=1:length(files)
        res = [res, load(['results' filesep files(i).name])];
    end
    
    % extract data
    rt_vals = [];
    for k=1:length(res)
        rt_vals = [rt_vals res(k).rT];
    end
    
    % Create table
    rt_vals = reshape(rt_vals,[3,6]);
    fm = [rt_vals(:,5) rt_vals(:,3) rt_vals(:,2) rt_vals(:,1) rt_vals(:,4) rt_vals(:,6)];
    T = table(fm);
    
    % Export to tex format
    if is_codeocean
        table2latex(T,'/results/logs/table5.tex')
    else
        table2latex(T,'table5.tex')
    end

end