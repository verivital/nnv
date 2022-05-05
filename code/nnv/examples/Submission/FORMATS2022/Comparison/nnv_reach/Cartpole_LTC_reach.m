%% Reachability analysis of the CTRNN Cartpole (12 dimensions)
% Idea is to compare to Gotube and other tools

% There is some dimension mismatch when trying to process the dynamics...
% Not sure what is happening... Try simply simulating the function using an
% ode solver

addpath('../benchmark_dynamics/');

reachstep = 0.001; % step size to compute reach sets
final_time = 0.001; % Time horizon
% MU = 1.3; % What is MU in Gotube?
Initial_radius = 1e-4; % Uncertainty in dynamics.
model = NonLinearODE(12,1,@CartPoleLTC, reachstep, final_time,eye(12));

% Change default options
% model.options.timeStep = 0.05;
% model.options.taylorTerms = 4;
model.options.zonotopeOrder = 100;
model.options.alg = 'lin';
model.options.tensorOrder = 5;
model.options.errorOrder3 = 10;

% Initial states
x0 = [0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
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
saveas(f,'../results/LTC_Cartpole_nnv.png');