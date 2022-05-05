%% Create MNIST results table
files_to_load = ["ffnn/ffnn_large_nnv","ffnn/ffnn_nnv","ffnn/ffnn_small_nnv",...
    "cnn/cnn_small_nnv","cnn/cnn_medium_nnv","cnn/cnn_tiny_nnv",...
    "ffnn/eval.mat","cnn/eval.mat"];

names = {'fnn_small';'fnn_mid';'fnn_large';'cnn_small';'cnn_mid';'cnn_large'};
% Run 1: inf 0.5
fnl = load(files_to_load{1}+"_inf_0.5.mat");
fnm = load(files_to_load{2}+"_inf_0.5.mat");
fns = load(files_to_load{3}+"_inf_0.5.mat");
cnl = load(files_to_load{4}+"_inf_0.5.mat");
cnm = load(files_to_load{5}+"_inf_0.5.mat");
cns = load(files_to_load{6}+"_inf_0.5.mat");
evalf = load(files_to_load{7});
evalc = load(files_to_load{8});
% Create tables
acc = [evalf.acc_small; evalf.acc_medium; evalf.acc_large; evalc.acc_tiny; evalc.acc_medium; evalc.acc_small];
robs = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT;cns.rob/cns.numT;cnm.rob/cnm.numT;cnl.rob/cnl.numT;];
timeA = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT;cns.timeT/cns.numT;cnm.timeT/cnm.numT;cnl.timeT/cnl.numT;];

% Run 2: inf 1
fnl = load(files_to_load{1}+"_inf_1.mat");
fnm = load(files_to_load{2}+"_inf_1.mat");
fns = load(files_to_load{3}+"_inf_1.mat");
cnl = load(files_to_load{4}+"_inf_1.mat");
cnm = load(files_to_load{5}+"_inf_1.mat");
cns = load(files_to_load{6}+"_inf_1.mat");

% Create columns
robs2 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT;cns.rob/cns.numT;cnm.rob/cnm.numT;cnl.rob/cnl.numT;];
timeA2 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT;cns.timeT/cns.numT;cnm.timeT/cnm.numT;cnl.timeT/cnl.numT;];

% Run 3: inf 2
fnl = load(files_to_load{1}+"_inf_2.mat");
fnm = load(files_to_load{2}+"_inf_2.mat");
fns = load(files_to_load{3}+"_inf_2.mat");
cnl = load(files_to_load{4}+"_inf_2.mat");
cnm = load(files_to_load{5}+"_inf_2.mat");
cns = load(files_to_load{6}+"_inf_2.mat");

% Create columns
robs3 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT;cns.rob/cns.numT;cnm.rob/cnm.numT;cnl.rob/cnl.numT;];
timeA3 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT;cns.timeT/cns.numT;cnm.timeT/cnm.numT;cnl.timeT/cnl.numT;];

% Run 4: random 1%
fnl = load(files_to_load{1}+"_random_2.55.mat");
fnm = load(files_to_load{2}+"_random_2.55.mat");
fns = load(files_to_load{3}+"_random_2.55.mat");
cnl = load(files_to_load{4}+"_random_2.55.mat");
cnm = load(files_to_load{5}+"_random_2.55.mat");
cns = load(files_to_load{6}+"_random_2.55.mat");

% Create columns
robs4 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT;cns.rob/cns.numT;cnm.rob/cnm.numT;cnl.rob/cnl.numT;];
timeA4 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT;cns.timeT/cns.numT;cnm.timeT/cnm.numT;cnl.timeT/cnl.numT;];

% Run 5: random 5%
fnl = load(files_to_load{1}+"_random_12.75.mat");
fnm = load(files_to_load{2}+"_random_12.75.mat");
fns = load(files_to_load{3}+"_random_12.75.mat");
cnl = load(files_to_load{4}+"_random_12.75.mat");
cnm = load(files_to_load{5}+"_random_12.75.mat");
cns = load(files_to_load{6}+"_random_12.75.mat");

% Create columns
robs5 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT;cns.rob/cns.numT;cnm.rob/cnm.numT;cnl.rob/cnl.numT;];
timeA5 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT;cns.timeT/cns.numT;cnm.timeT/cnm.numT;cnl.timeT/cnl.numT;];

% Run 6: random 10%
fnl = load(files_to_load{1}+"_random_25.5.mat");
fnm = load(files_to_load{2}+"_random_25.5.mat");
fns = load(files_to_load{3}+"_random_25.5.mat");
cnl = load(files_to_load{4}+"_random_25.5.mat");
cnm = load(files_to_load{5}+"_random_25.5.mat");
cns = load(files_to_load{6}+"_random_25.5.mat");

% Create columns
robs6 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT;cns.rob/cns.numT;cnm.rob/cnm.numT;cnl.rob/cnl.numT;];
timeA6 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT;cns.timeT/cns.numT;cnm.timeT/cnm.numT;cnl.timeT/cnl.numT;];

% names = ["NODEsmall";"NODEmid";"NODElarge";"CNODEsmall";"CNODEmid";"CNODElarge"];
% colmns = {'Accuracy (%)' 'Robust (10%)' 'Time (s)' 'Robust (10%)' 'Time (s)'}25.5; % Add another run of results?
T = table(acc,robs,timeA,robs2,timeA2,robs3,timeA3,robs4,timeA4,robs5,timeA5,robs6,timeA6);

% T.Properties.VariableNames = colmns;

table2latex(T,'mnist_res.tex')