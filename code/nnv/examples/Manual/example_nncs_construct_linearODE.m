
% /* An example of constructing a continuous linear plant model */
A = [0 1;-5 -2]; % system matrix
B = [0;3];       % control matrix
C = [0 1];       % output feedback matrix
D = [];          % output control matrix
Tc = 0.1;        % control period
Nr = 20; % number of reachability steps in one control period
sys = LinearODE(A, B, C, D, Tc, Nr); % plant object
