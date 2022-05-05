%% Reachability analysis of the 17-dimensional quadcopter
% Idea is to compare to Gotube and other tools

addpath('../benchmark_dynamics/');

reachstep = 1e-4; % step size to compute reach sets
final_time = 2; % Time horizon
% MU = 1.3; % What is MU in Gotube?
Initial_radius = 0.005; % Uncertainty in dynamics.
model = NonLinearODE(17,1,@Quadcopter, reachstep, final_time,eye(17));

% Change default options
% model.options.timeStep = 0.05;
% model.options.taylorTerms = 4;
% model.options.zonotopeOrder = 50;
model.options.alg = 'lin-adaptive';
% model.options.tensorOrder = 3;

% Initial states
x0 = [-0.995; -0.995; 9.005; -0.995; -0.995; -0.995; -0.995; -0.995; -0.995; 0; 0; 0; 1; 0; 0; 0; 0];
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
saveas(f,'../results/Quadcopter_nnv.png');