% if ~exist('../results','dir')
%     mkdir('../results')
% end
allTime = tic;
%% ACC
cd 'ACC';
run reachACC;
%% Airplane
cd ..;
cd 'Airplane'
run reach_airplane;
%% Benchmark 9 - Tora
cd ..;
cd 'Benchmark9-Tora';
run reach9;
%% Benchmark 10 - Unicycle
% cd ..;
% cd 'Benchmark10-Unicycle';
% run reach10.m;
%% Double Pendulum
% cd ..;
cd 'Double_Pendulum';
run reach_dp_less.m;
run reach_dp_more.m;
%% Single Pendulum
cd ..;
cd 'Single_Pendulum';
run reach_sp.m;
%% Tora Heterogeneous
cd ..;
cd 'Tora_Heterogeneous';
run reachTora_sigmoid.m;
run reachTora_reluTanh.m;
%% VCAS
cd ..;
cd VCAS;
run reachVCAS_middle19.m;
run reachVCAS_middle22.m;
run reachVCAS_middle25.m;
run reachVCAS_middle28.m;
run reachVCAS_worst19.m;
run reachVCAS_worst22.m;
run reachVCAS_worst25.m;
run reachVCAS_worst28.m;
%% Record overall time
allTime = toc(allTime);