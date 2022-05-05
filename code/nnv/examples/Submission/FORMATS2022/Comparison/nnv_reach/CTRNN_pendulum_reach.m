%% Reachability analysis of the CTRNN pendulum(10 dimensions)
% Idea is to compare to Gotube and other tools

addpath('../benchmark_dynamics/');

reachstep = 0.0005; % step size to compute reach sets
final_time = 0.05; % Time horizon
% MU = 1.3; % What is MU in Gotube?
Initial_radius = 1e-4; % Uncertainty in dynamics.
model = NonLinearODE(10,1,@PendulumwithCTRNN, reachstep, final_time,eye(10));

% Change default options
% model.options.timeStep = 0.05;
% model.options.taylorTerms = 4; % works
model.options.zonotopeOrder = 50; % works (100 - 3314s)
% model.options.alg = 'lin-adaptive';
% model.options.tensorOrder = 3;

% Initial states
x0 = zeros(10,1)+0.01;
lb = x0 - Initial_radius;
ub = x0 + Initial_radius;
init_set = Star(lb,ub);
input_set = Star(0,0); % No inputs, but need to define it

% Compute reachability analysis
t = tic;
R = model.stepReachStar(init_set,input_set);
toc(t);

% Plot reachable sets
f = figure;
Star.plotBoxes_2D_noFill(model.intermediate_reachSet,1,2,'b');
grid;
hold on;
% for sets = model.intermediate_reachSet
%     Zono.plots(sets.Z); % This is super slow, try saving using CORA
%     zonotopes and use their plotting routines (Tho this plots the actual reach set, not an overapproximation)
% end
xlabel('x1');
ylabel('x2');
saveas(f,'../results/CTRNN_Pendulum_nnv.png');