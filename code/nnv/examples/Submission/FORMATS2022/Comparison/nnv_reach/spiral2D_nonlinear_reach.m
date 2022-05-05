%% Reachability analysis of a spiral2D Neural Network ODE
% Function defined in a different file for CORA
C = eye(2); 
reachstep = 0.01; % step size to compute reach sets
final_time = 10.0; % Time horizon
Initial_radius = 0.01; % Uncertainty in dynamics (initial states)
% Contruct NeuralODE (odeblock)
model = NonLinearODE(2,1,@spiral_non,reachstep,final_time,C); % Nonlinear ODE plant 

% Change default options
% model.options.timeStep = 0.05;
% model.options.taylorTerms = 4;
% model.options.zonotopeOrder = 50;
% model.options.alg = 'lin-adaptive';
% model.options.tensorOrder = 3;

% Initial states
x0 = [2.0;0.0]; % This is like the initial input to the ODEblock (initial state)
lb = x0 - Initial_radius;
ub = x0 + Initial_radius;
init_set = Star(lb,ub);
input_set = Star(0,0); % No inputs, but need to define it

% Compute reachability analysis
t = tic;
R = model.stepReachStar(init_set,input_set);
time = toc(t);
Rall = model.intermediate_reachSet;
save('../results/spiral_nl_0.01.mat','Rall','time')

% Plot results
f = figure;
hold on;
Star.plotBoxes_2D_noFill(Rall,1,2,'b');
xlabel('x_1');
ylabel('x_2');
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 15; % Set font size of axis
% ax.YAxis.FontSize = 15;
% saveas(f,'nonlinearspiral_0.05.png')


%% Reachability scenario 2
% C = eye(2); 
% reachstep = 0.01; % step size to compute reach sets
% final_time = 10.0; % Time horizon
Initial_radius = 0.05; % Uncertainty in dynamics (initial states)
% Contruct NeuralODE (odeblock)
model = NonLinearODE(2,1,@spiral_non,reachstep,final_time,C); % Nonlinear ODE plant 

% Change default options
% model.options.timeStep = 0.05;
% model.options.taylorTerms = 4;
% model.options.zonotopeOrder = 50;
% model.options.alg = 'lin-adaptive';
% model.options.tensorOrder = 3;

% Initial states
% x0 = [10.0;0.0]; % This is like the initial input to the ODEblock (initial state)
lb = x0 - Initial_radius;
ub = x0 + Initial_radius;
init_set = Star(lb,ub);
input_set = Star(0,0); % No inputs, but need to define it

% Compute reachability analysis
t = tic;
R = model.stepReachStar(init_set,input_set);
time = toc(t);
Rall = model.intermediate_reachSet;
save('../results/spiral_nl_0.05.mat','Rall','time')

% Plot results
f = figure;
hold on;
Star.plotBoxes_2D_noFill(Rall,1,2,'b');
xlabel('x_1');
ylabel('x_2');
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 15; % Set font size of axis
% ax.YAxis.FontSize = 15;
% saveas(f,'nonlinearspiral_0.05.png')

%% Reachability scenario #3
% C = eye(2); 
% reachstep = 0.01; % step size to compute reach sets
% final_time = 10.0; % Time horizon
Initial_radius = 0.1; % Uncertainty in dynamics (initial states)
% Contruct NeuralODE (odeblock)
model = NonLinearODE(2,1,@spiral_non,reachstep,final_time,C); % Nonlinear ODE plant 

% Change default options
% model.options.timeStep = 0.05;
% model.options.taylorTerms = 4;
% model.options.zonotopeOrder = 50;
% model.options.alg = 'lin-adaptive';
% model.options.tensorOrder = 3;

% Initial states
% x0 = [2.0;0.0]; % This is like the initial input to the ODEblock (initial state)
lb = x0 - Initial_radius;
ub = x0 + Initial_radius;
init_set = Star(lb,ub);
input_set = Star(0,0); % No inputs, but need to define it

% Compute reachability analysis
t = tic;
R = model.stepReachStar(init_set,input_set);
time = toc(t);
Rall = model.intermediate_reachSet;
save('../results/spiral_nl_0.1.mat','Rall','time')

% Plot results
f = figure;
hold on;
Star.plotBoxes_2D_noFill(Rall,1,2,'b');
xlabel('x_1');
ylabel('x_2');
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 15; % Set font size of axis
% ax.YAxis.FontSize = 15;
% saveas(f,'nonlinearspiral_0.05.png')