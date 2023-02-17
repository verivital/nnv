
% /* An example of constructing an discrete linear NNCS object */
% /* construct a FFNN controller
L1 = LayerS([2; 1], [0.5; -1], 'poslin');
L2 = LayerS([1 -1], 0.2, 'purelin');
F  = FFNNS([L1 L2]);

% /* construct a plant model
A = [0 1;-5 -2]; % system matrix
B = [0;3];       % control matrix
C = [0 1];       % output feedback matrix
D = [];          % output control matrix
Ts = 0.1;        % sampling time
sys = DLinearODE(A, B, C, D, Ts); % plant object

% /* construct a linear NNCS object
ncs = DLinearNNCS(F, sys);