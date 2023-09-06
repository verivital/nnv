%% Reachability analysis of TORA
% Translational Oscillator with a Rotational Actuator (TORA)
% 
% This benchmark considers translational oscillations by a rotational actuator (TORA),
% where a cart is attached to a wall with a spring and is free to move on a friction-less surface.
% The cart has a weight attached to an arm inside it, which is free to
% rotate about an axis. 
% This benchark is part of the systems evaluate at the AINNCS ARCH friendly competition

% Load components and set reachability parameters
net = load_NN_from_mat('controllerTora.mat');
reachStep = 0.1; % can play with this
controlPeriod = 1; % fixed
plant = NonLinearODE(4,1,@dynamics9, reachStep, controlPeriod, eye(4));

% When do we want to finish the verification/simulation?
tFinal = 20;

% Initial set
lb = [0.6; -0.7; -0.4; 0.5];
ub = [0.7; -0.6; -0.3; 0.6];


%% Reachability analysis

% Create initial set
init_set = Star(lb,ub);

% Store all reachable sets
reachAll = init_set;

% Execute reachabilty analysis
reachOptions.reachMethod = 'approx-star';
alg = 'lin'; % plant reachability
% Begin reachability analysis (can do like this or using the NNCS class)
t = tic;
for i = 1:tFinal
    % Compute controller output set
    input_set = net.reach(init_set,reachOptions);
    input_set = input_set.getBox(); % overapprox the control action set due to CORA errors when too many constraints are in zonotope
    input_set = input_set.toStar();
    % Compute plant reachable set
    init_set = plant.stepReachStar(init_set, input_set, alg);
    reachAll = [reachAll init_set];
end
toc(t);

%% Visualize results
plant.get_interval_sets;

f = figure;
grid;
hold on;
rectangle('Position',[-2,-2,4,4],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
xlabel('x1');
ylabel('x2');

f1 = figure;
grid;
hold on;
rectangle('Position',[-2,-2,4,4],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
grid;
xlabel('x3');
ylabel('x4');

