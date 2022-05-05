%% Reachability analysis of the CTRNN Damped Forced Pendulum
% Idea is to compare to Gotube and other tools

addpath('../benchmark_dynamics/');

reachstep = 0.01; % step size to compute reach sets
final_time = 0.5; % Time horizon
% MU = 1.3; % What is MU in Gotube?
Initial_radius = 0.01; % Uncertainty in dynamics.
model = NonLinearODE(5,1,@CTRNN_DampedForcedPendulum, reachstep, final_time,eye(5));

% Change default options
% model.options.timeStep = 0.05;
% model.options.taylorTerms = 4;
% model.options.zonotopeOrder = 50;
% model.options.alg = 'lin-adaptive';
% model.options.tensorOrder = 2;

% Initial states
x0 = [0.21535, -0.58587, 0.8, 0.52323, 0.5]';
lb = x0 - Initial_radius;
ub = x0 + Initial_radius;
init_set = Star(lb,ub);
input_set = Star(0,0); % No inputs, but need to define it

% Compute reachability analysis
t = tic;
R = model.stepReachStar(init_set,input_set);
time = toc(t);
Rall = model.intermediate_reachSet;

save('../results/dfp_reach_short.mat','Rall','time')

% Plot reachable sets
f = figure;
Star.plotBoxes_2D_noFill(Rall,1,2,'b');
grid;
hold on;
% for sets = model.intermediate_reachSet
%     Zono.plots(sets.Z); % This is super slow, try saving using CORA
%     zonotopes and use their plotting routines (Tho this plots the actual reach set, not an overapproximation)
% end
xlabel('x1');
ylabel('x2');
saveas(f,'../results/CTRNN_DFP_nnv_short.png');