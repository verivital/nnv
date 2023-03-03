% Run all nnv reachability experiments
%% Set folder path for results
if ~exist('../nnvresults', 'dir')
       mkdir('../nnvresults')
end
addpath('benchmark_dynamics/');
% Turn off figure display
set(0,'DefaultFigureVisible','off')

%% Run benchmarks tool compare

% Spiral
run spiral2D_linear_reach.m

%% FPA
run CTRNN_FPA_short.m
run CTRNN_FPA_mid.m
run CTRNN_FPA_reach.m
