% This is a toy example in the paper
% Simulation-Equivalent Reachability of Large Linear Systems with Inputs
% Stanley Bak, CAV2017
A = [0 1; -1 0];
B = [1 0; 0 1];
lb = [-6; 0];
ub = [-5; 1];
Bi = Box(lb, ub);
X0 = Bi.toStar(); % input set: x[1] \in [-6, -5], x[2] \in [0, 1]
lb = [-0.5; -0.5];
ub = [0.5; 0.5];
Bu = Box(lb, ub);
U = Bu.toStar(); % control input set: -0.5 <= u[1], u[2] <= 0.5

h = pi/4; 
N = 2;

sys = LinearODE(A, B, [], []);
R1 = sys.simReach('direct', X0, U, h, N, []);
R2 = sys.simReach('ode45', X0, U, h, N, []);
R3 = sys.simReach('krylov', X0, U, h, N, 2);
figure;
Star.plots(R1);
figure;
Star.plots(R2);
figure;
Star.plots(R3);