clc;clear
load controller_Lcontainer_3in.mat; %load controller (mat format)
% Generate the controller object
n = length(W);
Layers = [];
for i=1:n
    L = Layer(W{i}, b{i}, layer_fcn(activation_fcns(i,:)));
    Layers = [Layers L];
end
Controller = FFNN(Layers); % feedforward neural network controller

Plant = NonLinearODE(9, 1, @LContainer);
Plant.set_timeStep(0.02); % time step for reachability analysis of the plant
Plant.set_tFinal(0.2); % Ts = 0.2, sampling time for control signal from neural network controller
output_mat = [0 0 1 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0 0]; % feedback 
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

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


[simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, []);

[P, reachTime] = ncs.reach('approx-star', init_set, input_ref, n_cores, N);

%fig = figure;
%Star.plotBoxes_2D_noFill(P, 1, 2, 'blue');
%saveas(fig, 'reachSet.pdf');
%save result.mat reachTime ncs;

%plot(simTrace(1, :), simTrace(2, :));