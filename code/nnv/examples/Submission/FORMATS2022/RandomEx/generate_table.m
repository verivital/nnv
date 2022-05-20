%% Generate time tables for paper submission
files = dir('results/*.mat'); % order alphabetically [L,M,S,XL,XS,XXL]
res = [];

for i=1:length(files)
    res = [res, load(['results' filesep files(i).name])];
end

% columns = [L,M,S,XL,XS,XXS];
% rows = [0.01;0.02;0.04];
rt_vals = [];
for k=1:length(res)
    rt_vals = [rt_vals res(k).rT];
end

rt_vals = reshape(rt_vals,[3,6]);
fm = [rt_vals(:,5) rt_vals(:,3) rt_vals(:,2) rt_vals(:,1) rt_vals(:,4) rt_vals(:,6)];
T = table(fm);

table2latex(T,'randEx_res.tex')