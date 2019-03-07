% This script performs rechability analysis for Emergency Braking System
% author: Dung Tran
% date: 6/3/2019

% normalization 
norm_mat = blkdiag([1/250 0 0; 0 3.6/120 0; 0 0 1/20], 3.6/120) ; % normalized matrix
norm_vec = [0; 0; 0; 0]; % normalized vector
norm_layer = LayerS(norm_mat, norm_vec, 'purelin');

% reinforcement learning controller
load controller.mat; 
rl_layer1 = LayerS(blkdiag(W{1, 1}, 1), [b{1, 1}'; 0], 'poslin'); % add one more neuron to represent v
rl_layer2 = LayerS(blkdiag(W{1, 2}, 1), [b{1, 2}'; 0], 'poslin'); % add one more neuron to represent v
rl_layer3 = LayerS(blkdiag(W{1, 3}, 1), [b{1, 3}';0], 'poslin'); % add one more neuron to represent v

% invert layer, first output of rl_layer3 is brake, second output is speed
% transformer layer requires first input is speed and second input is brake
invert_layer = LayerS([0 1; 1 0], [0;0], 'purelin');

% transformation 
load transform.mat;
tf_layer1 = LayerS([W{1, 1}; [0 1]], [b{1, 1}';0], 'poslin');
tf_layer2 = LayerS(blkdiag(W{1, 2},1), [b{1, 2}';0], 'poslin');
tf_layer3 = LayerS(blkdiag(W{1, 3},1), [b{1, 3}';0], 'purelin');


% control signal scale
scale_mat = [15*10 -15*120/3.6];
scale_vec = 0;
scale_layer = LayerS(scale_mat, scale_vec, 'purelin');


% composed controller

controller = FFNNS([norm_layer rl_layer1 rl_layer2 rl_layer3 invert_layer tf_layer1 tf_layer2 tf_layer3 scale_layer]);

% plant matrices
A = [1 -1/15 0; 0 1 0; 0 0 0];
B = [0;1/15;1];
C = [1 0 0;0 1 0; 0 0 1];
D = [0;0;0];
Ts = 0.1;
plant = DLinearODE(A, B, C, D, Ts);


% computation steps
% step 0: initial state of plant
% step 1: normalizing state
% step 2: compute brake output of rl_controller
% step 3: compose brake output and normalized speed
% step 4: compute transformer output
% step 5: scale control output
% step 6: compute reachable set of the plant with the new control input and
% initial state
% step 7: go back to step 1, ...


lb = [97; 25.2; -3.5];
ub = [97.5; 25.7; -3];
init_set = Star(lb, ub); % initial condition of the plant

N = 1; % number of control steps
numCores = 4;

X0 = init_set; % step 0: initial state of plant
S = X0;

%for i=1:N
speed = X0.affineMap([0 1 0], []);
X1 = Star.concatenateStars([X0 speed]);
L1 = norm_layer.reach(X1,'star',[]);
B1 = L1.getBox;
L2 = rl_layer1.reach(L1, 'star', []);
B2 = Star.get_hypercube_hull(L2);
L3 = rl_layer2.reach(L2, 'star', []);
B3 = Star.get_hypercube_hull(L3);
L4 = rl_layer3.reach(L3, 'star', []);
B4 = Star.get_hypercube_hull(L4);

L41 = invert_layer.reach(L4, 'star', []);
B41 = Star.get_hypercube_hull(L41);

L5 = tf_layer1.reach(L41, 'star', []);
B5 = Star.get_hypercube_hull(L5);

L6 = tf_layer2.reach(L5, 'star', []);
B6 = Star.get_hypercube_hull(L6);
L7 = tf_layer3.reach(L6, 'star', []);
B7 = Star.get_hypercube_hull(L7);
L8 = scale_layer.reach(L7, 'star', []);
B8 = Star.get_hypercube_hull(L8);

L9 = plant.stepReachStar(X0, B8.toStar);
B9 = Star.get_hypercube_hull(L9);

control = controller.reach(X1);
%control = Star.get_hypercube_hull(control);
%control = control.toStar;
%X0 = plant.stepReachStar(X0, control);

%end
