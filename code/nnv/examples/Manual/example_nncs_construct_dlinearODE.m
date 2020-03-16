
% /* An example of constructing a discrete linear plant model */
A = [0 1;-5 -2]; % system matrix
B = [0;3];       % control matrix
C = [0 1];       % output feedback matrix
D = [];          % output control matrix
Ts = 0.1;        % sampling time
sys = DLinearODE(A, B, C, D, Ts); % plant object
