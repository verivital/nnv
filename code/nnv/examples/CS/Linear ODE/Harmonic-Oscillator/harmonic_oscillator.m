% This is a toy example in the paper
% Simulation-Equivalent Reachability of Large Linear Systems with Inputs
% Stanley Bak, CAV2017
A = [0 1; -1 0];
B = [1 0; 0 1];

lb = [-6; 0];
ub = [-5; 1];
X0 = Star(lb,ub); 
lb = [-0.5; -0.5];
ub = [0.5; 0.5];
U = Star(lb,ub);

h = pi/4;
N = 3; % number of steps
contP = N*h; % control period

u = [0;0];
x0 = [-5.5;0.5];
t = [0:h:N];

sys = LinearODE(A, B, [], []);
[sy1,st1,sx1] = sys.simulate(repmat(u,[1,size(t,2)]),[0:h:N],x0);
R1 = sys.simReach('direct', X0, U, h, N, []);
R2 = sys.simReach('ode45', X0, U, h, N, []);
R3 = sys.simReach('krylov', X0, U, h, N, 2);
figure;
Star.plots(R1);
hold on;
plot(sx1(:,1),sx1(:,2),'--r');
figure;
Star.plots(R2);
hold on;
plot(sx1(:,1),sx1(:,2),'--r');
figure;
Star.plots(R3);
hold on;
plot(sx1(:,1),sx1(:,2),'--r');


%% Do the same using CORA linear models
lb = [-6; 0];
ub = [-5; 1];
X0 = Star(lb,ub); 
lb = [-0.5; -0.5];
ub = [0.5; 0.5];
U = Star(lb,ub);

t = [0:h:N];
C = eye(2);
D = [0 0;0 0];
sys2 = LinearODE_cora(A, B, C, D, h, contP);
R4 = sys2.stepReachStar(X0, U);
R4 = sys2.intermediate_reachSet;

figure;
Star.plots(R4);
hold on;
plot(sx1(:,1),sx1(:,2),'--r');
