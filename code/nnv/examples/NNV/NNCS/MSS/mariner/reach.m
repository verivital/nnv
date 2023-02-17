load controller_mariner_3in.mat; %load controller (mat format)
% Generate the controller object
n = length(W);
Layers = [];
for i=1:n
    L = LayerS(W{i}, b{i}, layer_fcn(activation_fcns(i,:)));
    Layers = [Layers L];
end
Controller = FFNNS(Layers); % feedforward neural network controller

reachStep = 0.02;
controlPeriod = 1;
output_mat = [0 0 1 0 0 0 0;
    0 0 0 0 0 1 0]; % feedback
Plant = NonLinearODE(7, 1, @Mariner, reachStep, controlPeriod, output_mat);

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial condition of the Plant
lb = [0; 0; 0; 0; 0; 0; 0];
ub = [0.0; 0.0; 0.0; 0.0; 0.0; 0.01; 0.0];
init_set = Star(lb, ub);

% reference input for neural network controller
lb1 = [5*pi/180];
ub1 = [5*pi/180];
input_ref = Star(lb1, ub1);

N = 10; % number of control steps   

n_cores = 4; % number of cores


[simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, lb1);

reachPRM.ref_input = input_ref;
reachPRM.init_set = init_set;
reachPRM.numSteps = N;
reachPRM.numCores = n_cores;
reachPRM.reachMethod = 'approx-star';
[P, reachTime] = ncs.reach(reachPRM);

fig = figure;
Star.plotBoxes_2D_noFill(P, 6, 6, 'blue');
%saveas(fig, 'reachSet.pdf');
%save result.mat reachTime ncs;

plot(simTrace(2, :), simTrace(1, :));