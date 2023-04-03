% function table_other();

% Process results of comparison except for acas

% Load data
cd oval21;
oval = load("results_oval21_approx_all.mat");
cd ..;

cd rl_benchmarks;
rl = load("results_rl.mat");
cd ..;

cd  tllverify;
tll = load("results_tllverify.mat");
cd ..;

% OVAL (no supported by MATLAB)
nnv_oval_time  = sum(oval.res(:,2))/30;
nnv_oval_sat   = sum(oval.res(:,1) == 0);
nnv_oval_unsat = sum(oval.res(:,1) == 1);

% RL benchmarks
mat_rl_time  = sum(rl.res(:,2))/50;
mat_rl_sat   = sum(rl.res(:,1) == 0);
mat_rl_unsat = sum(rl.res(:,1) == 1);
nnv_rl_time  = sum(rl.res(:,4))/50;
nnv_rl_sat   = sum(rl.res(:,3) == 0);
nnv_rl_unsat = sum(rl.res(:,3) == 1);

% tllverify
mat_tll_time  = sum(tll.res(:,2))/24;
mat_tll_sat   = sum(tll.res(:,1) == 0);
mat_tll_unsat = sum(tll.res(:,1) == 1);
nnv_tll_time  = sum(tll.res(:,4))/24;
nnv_tll_sat   = sum(tll.res(:,3) == 0);
nnv_tll_unsat = sum(tll.res(:,3) == 1);

% end