% Reachability analysis of a Neural Network ODE

%% Define layers and neural ODE
controlPeriod = 10; % total seconds
reachStep = 0.01; % 1 second
C = eye(2); % Want to get both of the outputs from NeuralODE
% Load parameters
load('odeffnn_spiral.mat');
% Contruct NeuralODE
% ODEBlock only linear layers
% Convert in form of a linear ODE model
states = 2;
outputs = 2;
inputs = 1;
w1 = Wb{3};
b1 = Wb{4}';
w2 = Wb{5};
b2 = Wb{6}';
Aout = w2*w1;
Bout = b2' + b1'*w2';
Cout = eye(states);
D = zeros(outputs,1);
numSteps = controlPeriod/reachStep;
odeblock = LinearODE(Aout,Bout',Cout,D,controlPeriod,numSteps);
odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
neuralode = NeuralODE({odelayer});

x0 = [2.0;0.0]; % This is like the initial input to the ODEblock (initial state)


%% Reachability run #1
% Setup
Initial_radius = 0.01; % Uncertainty in dynamics (initial states)
lb = x0 - Initial_radius;
ub = x0 + Initial_radius;
R0 = Star(lb,ub);

t = tic;
Rall = neuralode.reach(R0); % Reachability
time = toc(t);
save('../results/spiral_0.01.mat','Rall','time')

% Plot results
f = figure;
hold on;
Star.plotBoxes_2D_noFill(Rall,1,2,'b');
xlabel('x_1');
ylabel('x_2');

%% Reachability run #1
% Setup
Initial_radius = 0.05; % Uncertainty in dynamics (initial states)
lb = x0 - Initial_radius;
ub = x0 + Initial_radius;
R0 = Star(lb,ub);

t = tic;
Rall = neuralode.reach(R0); % Reachability
time = toc(t);
save('../results/spiral_0.05.mat','Rall','time')

% Plot results
f = figure;
hold on;
Star.plotBoxes_2D_noFill(Rall,1,2,'b');
xlabel('x_1');
ylabel('x_2');

%% Reachability run #3
% Setup
Initial_radius = 0.1; % Uncertainty in dynamics (initial states)
lb = x0 - Initial_radius;
ub = x0 + Initial_radius;
R0 = Star(lb,ub);

t = tic;
Rall = neuralode.reach(R0); % Reachability
time = toc(t);
save('../results/spiral_0.1.mat','Rall','time')

% Plot results
f = figure;
hold on;
Star.plotBoxes_2D_noFill(Rall,1,2,'b');
xlabel('x_1');
ylabel('x_2');