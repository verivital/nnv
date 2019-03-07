% This script performs rechability analysis for Emergency Braking System
% author: Dung Tran
% date: 6/3/2019

% normalization 
norm_mat = [1/250 0 0; 0 3.6/120 0; 0 0 1/20]; % normalized matrix

% reinforcement learning controller
load controller.mat; 
rl_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
rl_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
rl_layer3 = LayerS(W{1, 3}, b{1, 3}', 'poslin');
rl_controller = FFNNS([rl_layer1 rl_layer2 rl_layer3]); 

% transformation 
load transform.mat;
tf_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
tf_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
tf_layer3 = LayerS(W{1, 3}, b{1, 3}', 'purelin');
transformer = FFNNS([tf_layer1 tf_layer2 tf_layer3]);

% control signal scale
scale_mat = [-15*120/3.6 15*10];

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
ub = [97.5; 25.5; -3.0];
init_set = Star(lb, ub); % initial condition of the plant

N = 3; % number of control steps
numCores = 4;

X0 = init_set; % step 0: initial state of plant
S = X0;

speed_brakes = [];
tf_outs = [];
for i=1:N
    % step 1: normalizing state
    norm_X = X0.affineMap(norm_mat, []); % normalized state
    % step 2: compute brake output of rl_controller
    brake = rl_controller.reach(norm_X);
    brake = Star.get_hypercube_hull(brake);
    % step 3: compose brake output and speed
    [speed_lb, speed_ub] = norm_X.getRange(2);
    speed_brake = Star([speed_lb; brake.lb], [speed_ub; brake.ub]);
    speed_brakes = [speed_brakes speed_brake.getBox];
    % step 4: compute transformer output
    tf_out = transformer.reach(speed_brake);
    tf_out = Star.get_hypercube_hull(tf_out);
    tf_outs = [tf_outs tf_out];
    % step 5: scale control output
    control = Star([speed_lb; tf_out.lb], [speed_ub; tf_out.ub]); 
    control = control.affineMap(scale_mat, []);   
    % step 6: compute reachable set of the plant with the new control input and
    % initial state
    X0 = plant.stepReachStar(X0, control);  
    S = [S X0];
end

figure;
Star.plotBoxes_2D_noFill(S, 1, 2, 'b');
xlabel('distance');
ylabel('speed');

