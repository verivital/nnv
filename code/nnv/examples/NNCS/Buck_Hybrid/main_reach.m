%% Reachability analysis of the buck converter using NNV
% Load hybrid automata NNCS
plant = HybridA(5,1,buck_v2,2);
% Create random controller
W{1} = single(rand(64,5));
W{2} = single(rand(64,64));
W{3} = single(rand(1,64));
b{1} = single(rand(1,64));
b{2} = single(rand(1,64));
b{3} = single(rand(1,1));
layer1 = LayerS(W{1},b{1}','poslin'); % hidden layer #1 (ReLU)
layer2 = LayerS(W{2},b{2}','poslin'); % hidden layer #2 (ReLU)
layer3 = LayerS(W{3},b{3}','satlins'); % output layer (satlins)
Layers = [layer1 layer2 layer3];
controller = FFNNS(Layers); % neural network controller
feedbackMap = [0];

steps = 5; % Number of control steps
Ts = 1.6667e-05 / 10; % 1/10 switching frequency
plant.set_tFinal(steps*Ts);
plant.set_timeStep(Ts);

nncs = HybridANNCS(controller,plant,feedbackMap);

%% Setup reachability parameters

% Execute reachability analysis star set NNV
inp_set = Star;
lb = [-0.09; 0; 0; 0.75; -0.09];
ub = [0.11; 0; 0; 0.75; 0.11];
init_set = Star(lb,ub);

reachPRM.reachMethod = 'approx-star';
reachPRM.init_set = init_set;
reachPRM.numCores = 1;
reachPRM.ref_input = [];
reachPRM.numSteps = steps;

R = nncs.reach(reachPRM); % Seems to work at the moment

%% Visualize results

Sall = nncs.plant.intermediate_reachSet;

figure;
Star.plotBoxes_2D_noFill(Sall,3,1,'b');
hold on;
xlabel('x_3');
ylabel('x_1');

figure;
Star.plotBoxes_2D_noFill(Sall,3,5,'b');
xlabel('x_3');
ylabel('x_5');

figure;
Star.plotBoxes_2D_noFill(Sall,1,5,'b');
xlabel('x_1');
ylabel('x_5');

figure;
Star.plotBoxes_2D_noFill(Sall,3,2,'b');
xlabel('x_3');
ylabel('x_2');

figure;
Star.plotBoxes_2D_noFill(Sall,4,5,'b');
xlabel('x_4');
ylabel('x_5');
