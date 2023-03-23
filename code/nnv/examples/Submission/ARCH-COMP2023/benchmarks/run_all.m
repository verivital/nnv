% Supress warnings
w = warning ('off','all');

% ACC 
cd ACC;
reach(); % verified (~ 15 seconds)
cd ..;

% Airplane
cd Airplane;
run reach.m; % Falsified
cd ..;

% Double Pendulum (final)
cd Double_Pendulum;
run reach_more.m; % falsified
run reach_less.m; % falsified
cd ..;

% Single Pendulum (final)
cd Single_Pendulum;
run reach.m; % verified
cd ..;

% VCAS
cd VCAS;
run run_vcas.m; % verified
cd ..;

% Tora
cd Tora_Heterogeneous;
reachTora_reluTanh(); % verified with input partition (~ 38 mins)
reachTora_sigmoid(); % verified with input partition (~ 2.7 hours, let's try a different partition to speed it up)
cd ..;
% cd Benchmark9-Tora;
% reach(); % unknown
% cd ..;

% Quad
% Unknown

% Attitude
% Unknown 

% Unicycle
% unknown -> overapproximation

% Docking 
% unknown -> overapproximation