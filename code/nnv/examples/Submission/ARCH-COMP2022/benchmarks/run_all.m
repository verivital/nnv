% ACC (final)
cd ACC;
run reach.m; % verified
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
run reachTora_reluTanh.m;
% sigmoid -> running
% relu+tanh -> verified
% Tora relu -> unknown 
% cd ..;

% Quad
% Unknown

% Attitude
% Unknown 

% Unicycle
% unknown -> overapproximation

% Docking 
% unknown -> overapproximation