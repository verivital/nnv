%% Create a random and example of a Neural ODE
rng(0);

%% Create random weight variables
% Set as global variables the parameters used in the nonlinear dynamics
states = 3;
input = 2;
output = 2;
h = 5;
global w21 b21 w22 b22 w4 b4;
w1 = randn(input,states)'; b1 = randn(states,1);     % layer 1
w21 = randn(h,states); b21 = randn(h,1); % layer 2.1
w22 = randn(states,h); b22 = randn(states,1);  % layer 2.2
w3 = randn(states,states)'; b3 = randn(states,1);   % layer 3
w4 = randn(states,states); b4 = randn(states,1);  % layer 4
w5 = randn(states,output)'; b5 = randn(output,1);   % layer 5

%% Create NeuralODE 
layer1 = LayerS(w1,b1,'tansig'); % tanh
layer3 = LayerS(w3,b3,'tansig'); % tanh
layer5 = LayerS(w5,b5,'tansig'); % tanh

cP1 = 0.1; % tf simulation
reachStep1 = 0.005; % reach step
C1 = eye(states);
odeblock1 = NonLinearODE(states,1,@dyn1,reachStep1,cP1,C1); % Nonlinear ODE
odeblock1.options.tensorOrder = 2;
layer2 = ODEblockLayer(odeblock1,cP1,reachStep1,false);
cP2 = 1; % tf simulation
reachStep2 = 0.01; % reach step
C2 = eye(states);
odeblock2 = NonLinearODE(states,1,@dyn2,reachStep2,cP2,C2);
odeblock2.options.tensorOrder = 2;
layer4 = ODEblockLayer(odeblock2,cP1,reachStep1,true);
neuralLayers = {layer1,layer2,layer3,layer4,layer5};
neuralode = NeuralODE(neuralLayers);
% Input
x0 = rand(input,1); 
nname = "S";

%% Reachability run #1 
unc = 0.01;
reach_random(neuralode,x0,nname,unc);
disp('Completed');

%% Reachability run #2 
% Reset the neuralode
odeblock1 = NonLinearODE(states,1,@dyn1,reachStep1,cP1,C1); % Nonlinear ODE
odeblock1.options.tensorOrder = 2;
layer2 = ODEblockLayer(odeblock1,cP1,reachStep1,false);
odeblock2 = NonLinearODE(states,1,@dyn2,reachStep2,cP2,C2);
odeblock2.options.tensorOrder = 2;
layer4 = ODEblockLayer(odeblock2,cP1,reachStep1,true);
neuralLayers = {layer1,layer2,layer3,layer4,layer5};
neuralode = NeuralODE(neuralLayers);
% Run
unc = 0.02;
reach_random(neuralode,x0,nname,unc);
disp('Completed');

%% Reachability run #3 
odeblock1 = NonLinearODE(states,1,@dyn1,reachStep1,cP1,C1); % Nonlinear ODE
odeblock1.options.tensorOrder = 2;
layer2 = ODEblockLayer(odeblock1,cP1,reachStep1,false);
odeblock2 = NonLinearODE(states,1,@dyn2,reachStep2,cP2,C2);
odeblock2.options.tensorOrder = 2;
layer4 = ODEblockLayer(odeblock2,cP1,reachStep1,true);
neuralLayers = {layer1,layer2,layer3,layer4,layer5};
neuralode = NeuralODE(neuralLayers);
% Run
unc = 0.04;
reach_random(neuralode,x0,nname,unc);
disp('Completed');

%% Dynamical functions
% Dynamics of first ODElayer
function dx = dyn1(x,t)
    global w21 b21 w22 b22;
    dx1 = tanh(w21*x+b21);
    dx = tanh(w22*dx1+b22);
end

function dx = dyn2(x,t)
    global w4 b4;
    dx = tanh(w4*x+b4);
end
