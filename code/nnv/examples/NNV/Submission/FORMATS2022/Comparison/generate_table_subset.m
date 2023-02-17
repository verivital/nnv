%% Partially reproduce table 3
% Get times from all the tools, then generate able
%% GoTube results 
% Time of these files are saved in a json file under the logged directory
% These files are simply numbered in the order they are produced, so
% assuming that the last file we run was the run_formats_subset one, we get
% the last 3 json files from this repo
gotubeLpath = 'GoTube/logged/';
if is_codeocean
    gotubeLpath = ['/results/logs/tuberesults/logged/'];
end
tubefiles = dir(gotubeLpath);
% Order of execution is
% 1) fpa3
% 2) fpa2
% 3) fpa1
% 4) spiralL1
% 5) spiralL2
% 6) spiralL3
% l3_gt = get_time_gotube([gotubeLpath, tubefiles(end).name]);
% l2_gt = get_time_gotube([gotubeLpath, tubefiles(end-1).name]);
% l1_gt = get_time_gotube([gotubeLpath, tubefiles(end-2).name]);
% fpa1_gt = get_time_gotube([gotubeLpath, tubefiles(end-3).name]);
% fpa2_gt = get_time_gotube([gotubeLpath, tubefiles(end-4).name]);
% fpa3_gt = get_time_gotube([gotubeLpath, tubefiles(end-5).name]);
fpa1_gt = get_time_gotube([gotubeLpath, tubefiles(end).name]);
fpa2_gt = get_time_gotube([gotubeLpath, tubefiles(end-1).name]);
fpa3_gt = get_time_gotube([gotubeLpath, tubefiles(end-2).name]);

%% Flowstar
flowLpath = 'flowstar/results/';
if is_codeocean
    flowLpath = ['/results/logs/flowresults/'];
end
l1_f = get_time_flowstar([flowLpath, 'spiral_L1.txt']);
l2_f = get_time_flowstar([flowLpath, 'spiral_L2.txt']);
l3_f = get_time_flowstar([flowLpath, 'spiral_L3.txt']);
% fpa1_f = get_time_flowstar([flowLpath, 'fpa_time_short.txt']);
% fpa2_f = get_time_flowstar([flowLpath, 'fpa_time_mid.txt']);
% fpa3_f = get_time_flowstar([flowLpath, 'fpa_time.txt']);

%% Juliareach
juliaLpath = 'juliareach/results/';
if is_codeocean
    juliaLpath = ['/results/logs/juliaresults/'];
end
spiral = load([juliaLpath, 'spiral_time.mat']);
l1_jl = spiral.time_4;
l2_jl = spiral.time_5;
l3_jl = spiral.time_6;
fpa_jl = load([juliaLpath, 'fpa_time.mat']);
fpa1_jl = fpa_jl.time_short;
fpa2_jl = fpa_jl.time_mid;
fpa3_jl = fpa_jl.time_2;

%% NNV
nnvLpath = 'nnvresults/';
% Doing this locally, then moving tables and figures to correct folder in
% codeOcean
l1_nnv = load([nnvLpath, 'spiral_0.01.mat']);
l1_nnv = l1_nnv.time;
l2_nnv = load([nnvLpath, 'spiral_0.05.mat']);
l2_nnv = l2_nnv.time;
l3_nnv = load([nnvLpath, 'spiral_0.1.mat']);
l3_nnv = l3_nnv.time;
fpa1_nnv = load([nnvLpath, 'fpa_reach_short.mat']);
fpa1_nnv = fpa1_nnv.time;
fpa2_nnv = load([nnvLpath, 'fpa_reach_mid.mat']);
fpa2_nnv = fpa2_nnv.time;
fpa3_nnv = load([nnvLpath, 'fpa_reach.mat']);
fpa3_nnv = fpa3_nnv.time;

%% Create table
flowstar = [string(l1_f); string(l2_f); string(l3_f); " -- "; " --"; "--"];
gotube = [" -- "; " --"; "--"; string(fpa1_gt); string(fpa2_gt); string(fpa3_gt)];
juliareach = [l1_jl; l2_jl; l3_jl; fpa1_jl; fpa2_jl; fpa3_jl];
nnvode = [l1_nnv; l2_nnv; l3_nnv; fpa1_nnv; fpa2_nnv; fpa3_nnv];
% generate table
T = table(flowstar, gotube, juliareach, nnvode)
if is_codeocean
    table2latex(T,'/results/logs/table3.tex');
else
    table2latex(T,'table3.tex');
end