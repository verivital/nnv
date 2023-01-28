function tool_comparison_table()

    %% Reproduce table 3
    % Get times from all the tools, then generate table
   
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
    % 4) spiralL1 % this fails, no files generated
    % 5) spiralL2 % this fails, no files generated
    % 6) spiralL3 % this fails, no files generated
    % 7) spiralNL1
    % 8) spiralNL2
    % 9) spiralNL3
    % 10) cartpole3
    % 11) cartpole2
    % 12) cartpole1 
    
    % Cartpole
    ct1_gt = get_time_gotube([gotubeLpath, tubefiles(end).name]);
    ct2_gt = get_time_gotube([gotubeLpath, tubefiles(end-1).name]);
    ct3_gt = get_time_gotube([gotubeLpath, tubefiles(end-2).name]);
    % Nonlinear Spiral
    nl3_gt = get_time_gotube([gotubeLpath, tubefiles(end-3).name]);
    nl2_gt = get_time_gotube([gotubeLpath, tubefiles(end-4).name]);
    nl1_gt = get_time_gotube([gotubeLpath, tubefiles(end-5).name]);
    % FPA
    fpa1_gt = get_time_gotube([gotubeLpath, tubefiles(end-6).name]);
    fpa2_gt = get_time_gotube([gotubeLpath, tubefiles(end-7).name]);
    fpa3_gt = get_time_gotube([gotubeLpath, tubefiles(end-8).name]);
    
    %% Flowstar

    % Load flowstar results
    flowLpath = 'flowstar/results/';
    if is_codeocean
        flowLpath = ['/results/logs/flowresults/'];
    end
    % Linear Spiral
    l1_f = get_time_flowstar([flowLpath, 'spiral_L1.txt']);
    l2_f = get_time_flowstar([flowLpath, 'spiral_L2.txt']);
    l3_f = get_time_flowstar([flowLpath, 'spiral_L3.txt']);
    % Nonlinear spiral
    % nl1_f = get_time_flowstar([flowLpath, 'spiral_NL1.txt']); % This and following 9 files are missing as flow* timed out
    % nl2_f = get_time_flowstar([flowLpath, 'spiral_NL2.txt']);
    % nl3_f = get_time_flowstar([flowLpath, 'spiral_NL3.txt']);
    % FPA
    % fpa1_f = get_time_flowstar([flowLpath, 'fpa_time_short.txt']);
    % fpa2_f = get_time_flowstar([flowLpath, 'fpa_time_mid.txt']);
    % fpa3_f = get_time_flowstar([flowLpath, 'fpa_time.txt']);
    % Cartpole
    % ct1_f = get_time_flowstar([flowLpath, 'cartpole_time_short.txt']);
    % ct2_f = get_time_flowstar([flowLpath, 'cartpole_time_mid.txt']);
    % ct3_f = get_time_flowstar([flowLpath, 'cartpole_time.txt']);
    

    %% Juliareach
    
    % Load juliareach results
    juliaLpath = 'juliareach/results/';
    if is_codeocean
        juliaLpath = ['/results/logs/juliaresults/'];
    end
    spiral = load([juliaLpath, 'spiral_time.mat']);
    % Nonlinear Spiral
    nl1_jl = spiral.time_1;
    nl2_jl = spiral.time_2;
    nl3_jl = spiral_time_3;
    % Linear Spiral
    l1_jl = spiral.time_4;
    l2_jl = spiral.time_5;
    l3_jl = spiral.time_6;
    % FPA
    fpa_jl = load([juliaLpath, 'fpa_time.mat']);
    fpa1_jl = fpa_jl.time_short;
    fpa2_jl = fpa_jl.time_mid;
    fpa3_jl = fpa_jl.time_2;
    % Cartpole
    ct_jl = load([juliaLpath, 'Cartpole_time.mat']);
    ct1_jl = ct_jl.time_short;
    ct2_jl = ct_jl.time_mid;
    ct3_jl = ct_jl.time_2;
    

    %% NNV

    % Doing this locally, then moving tables and figures to correct folder in CodeOcean/Docker
    nnvLpath = 'nnvresults/';

    % Linear Spiral
    l1_nnv = load([nnvLpath, 'spiral_0.01.mat']);
    l1_nnv = l1_nnv.time;
    l2_nnv = load([nnvLpath, 'spiral_0.05.mat']);
    l2_nnv = l2_nnv.time;
    l3_nnv = load([nnvLpath, 'spiral_0.1.mat']);
    l3_nnv = l3_nnv.time;
    % Nonlinear spiral
    nl1_nnv = load([nnvLpath, 'spiral_nl_0.01.mat']);
    nl1_nnv = nl1_nnv.time;
    nl2_nnv = load([nnvLpath, 'spiral_nl_0.05.mat']);
    nl2_nnv = nl2_nnv.time;
    nl3_nnv = load([nnvLpath, 'spiral_nl_0.1.mat']);
    nl3_nnv = nl3_nnv.time;
    % FPA
    fpa1_nnv = load([nnvLpath, 'fpa_reach_short.mat']);
    fpa1_nnv = fpa1_nnv.time;
    fpa2_nnv = load([nnvLpath, 'fpa_reach_mid.mat']);
    fpa2_nnv = fpa2_nnv.time;
    fpa3_nnv = load([nnvLpath, 'fpa_reach.mat']);
    fpa3_nnv = fpa3_nnv.time;
    % Cartpole
    ct1_nnv = load([nnvLpath, 'cartpole_reach_short.mat']);
    ct1_nnv = ct1_nnv.time;
    ct2_nnv = load([nnvLpath, 'cartpole_reach_mid.mat']);
    ct2_nnv = ct2_nnv.time;
    ct3_nnv = load([nnvLpath, 'cartpole_reach.mat']);
    ct3_nnv = ct3_nnv.time;
    
    %% Create table
    
    % Rows are: spiralL (1-3), spiralNL (4-6), FPA (7-9), Cartpole (10-12)
    flowstar = [string(l1_f); string(l2_f); string(l3_f); " -- "; " --"; "--"; " -- "; " --"; "--"; " -- "; " --"; "--"];
    gotube = [" -- "; " --"; "--"; string(nl1_gt); string(nl2_gt); string(nl3_gt); string(fpa1_gt); string(fpa2_gt); string(fpa3_gt); string(ct1_gt); string(ct2_gt); string(ct3_gt)];
    juliareach = [l1_jl; l2_jl; l3_jl; nl1_jl; nl2_jl; nl3_jl; fpa1_jl; fpa2_jl; fpa3_jl; ct1_jl; ct2_jl; ct3_jl;];
    nnvode = [l1_nnv; l2_nnv; l3_nnv; nl1_nnv; nl2_nnv; nl3_nnv; fpa1_nnv; fpa2_nnv; fpa3_nnv; ct1_nnv; ct2_nnv; ct3_nnv;];
    
    % generate table and save it
    T = table(flowstar, gotube, juliareach, nnvode);
    if is_codeocean
        table2latex(T,'/results/logs/table3.tex');
    else
        table2latex(T,'table3.tex');
    end