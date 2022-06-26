%% Create MNIST results table
files_to_load = ["ffnn/ffnn_large_nnv","ffnn/ffnn_nnv","ffnn/ffnn_small_nnv","ffnn/eval.mat"];

names = {'fnn_small';'fnn_mid';'fnn_large'};
% Run 1: inf 0.5
fnl = load(files_to_load{1}+"_inf_0.5.mat");
fnm = load(files_to_load{2}+"_inf_0.5.mat");
fns = load(files_to_load{3}+"_inf_0.5.mat");

evalf = load(files_to_load{4});

% Create tables
acc = [evalf.acc_small; evalf.acc_medium; evalf.acc_large];
robs = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT];
timeA = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT];

% Run 2: inf 1
fnl = load(files_to_load{1}+"_inf_1.mat");
fnm = load(files_to_load{2}+"_inf_1.mat");
fns = load(files_to_load{3}+"_inf_1.mat");

% Create columns
robs2 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT];
timeA2 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT];

% Run 3: inf 2
fnl = load(files_to_load{1}+"_inf_2.mat");
fnm = load(files_to_load{2}+"_inf_2.mat");
fns = load(files_to_load{3}+"_inf_2.mat");

% Create columns
robs3 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT];
timeA3 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT];

% Run 4: random 1%
fnl = load(files_to_load{1}+"_random_2.55.mat");
fnm = load(files_to_load{2}+"_random_2.55.mat");
fns = load(files_to_load{3}+"_random_2.55.mat");

% Create columns
robs4 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT];
timeA4 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT];

% Run 5: random 5%
fnl = load(files_to_load{1}+"_random_12.75.mat");
fnm = load(files_to_load{2}+"_random_12.75.mat");
fns = load(files_to_load{3}+"_random_12.75.mat");

% Create columns
robs5 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT];
timeA5 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT];

% Run 6: random 10%
fnl = load(files_to_load{1}+"_random_25.5.mat");
fnm = load(files_to_load{2}+"_random_25.5.mat");
fns = load(files_to_load{3}+"_random_25.5.mat");

% Create columns
robs6 = [fns.rob/fns.numT;fnm.rob/fnm.numT;fnl.rob/fnl.numT];
timeA6 = [fns.timeT/fns.numT;fnm.timeT/fnm.numT;fnl.timeT/fnl.numT];

% names = ["NODEsmall";"NODEmid";"NODElarge";"CNODEsmall";"CNODEmid";"CNODElarge"];
% colmns = {'Accuracy (%)' 'Robust (10%)' 'Time (s)' 'Robust (10%)' 'Time (s)'}25.5; % Add another run of results?
T = table(acc,robs,timeA,robs2,timeA2,robs3,timeA3,robs4,timeA4,robs5,timeA5,robs6,timeA6);

% T.Properties.VariableNames = colmns;

if is_codeocean
    table2latex(T,'/results/logs/table4.tex')
else
    table2latex(T,'table4.tex')
end