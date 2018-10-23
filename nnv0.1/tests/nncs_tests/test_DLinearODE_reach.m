A = [0 1;-5 -2];
B = [0;3];
C = [0 1];
D = 0;
Ts = 0.1;
sys = LinearODE(A, B, C, D);
sysd = sys.c2d(Ts); % convert from continuous to discrete

% set of initial condition
%      -1 <= x0[1] <= 1; -1 <= x0[2] <= 1

lb1 = [-1 ;-1];
ub1 = [1; 1]; 
B = Box(lb1, ub1);
I = B.toStar();

% set of control input: -0.5 <= u <= 1
lb1 = [-1];
ub1 = [0.5];
B = Box(lb1, ub1);
U = B.toStar();

fig = figure;
U.plot; % plot control input set
fig = figure;
I.plot; % plot initial condition


R = sysd.stepReachStar(I, U);
fig = figure;
R.plot; % plot one step reachable set

