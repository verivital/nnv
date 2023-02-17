% clc;clear
load controller_Lcontainer_3in.mat; %load controller (mat format)
% Generate the controller object
n = length(W);
Layers = [];
for i=1:n
    L = LayerS(W{i}, b{i}, layer_fcn(activation_fcns(i,:)));
    Layers = [Layers L];
end
Controller = FFNNS(Layers); % feedforward neural network controller

reachStep = 0.02;
controlPeriod = 0.2;
output_mat = [0 0 1 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0 0]; % feedback 
Plant = NonLinearODE(9, 1, @LContainer, reachStep, controlPeriod, output_mat);

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial condition of the Plant
lb = [0; 0; 0; 0; 0; 0; 0; 0; 0];
ub = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1];
init_set = Star(lb, ub);

% reference input for neural network controller
lb1 = [5*pi/180];
ub1 = [5*pi/180];
input_ref = Star(lb1, ub1);

N = 2; % number of control steps   

n_cores = 4; % number of cores

% [simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, []);

reachPRM.init_set = init_set;
reachPRM.ref_input = input_ref;
reachPRM.numCores = n_cores;
reachPRM.numSteps = N;
reachPRM.reachMethod = 'approx-star';
% [P, reachTime] = ncs.reach(reachPRM); % Error when creating the hessian
% tensors in CORA

%fig = figure;
%Star.plotBoxes_2D_noFill(P, 1, 2, 'blue');
%saveas(fig, 'reachSet.pdf');
%save result.mat reachTime ncs;

%plot(simTrace(1, :), simTrace(2, :));