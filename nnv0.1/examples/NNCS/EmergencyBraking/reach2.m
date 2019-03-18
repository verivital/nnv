% This script performs rechability analysis for Emergency Braking System
% with the RL controller that has saturation activation function
% author: Dung Tran
% date: 6/3/2019

% normalization 
norm_mat = [1/250 0 0; 0 3.6/120 0; 0 0 1/20]; % normalized matrix

% reinforcement learning controller
load ControllerSat.mat; 
rl_layer1 = LayerS(controller.W{1, 1}, controller.b{1, 1}', 'poslin');
rl_layer2 = LayerS(controller.W{1, 2}, controller.b{1, 2}', 'poslin');
rl_layer3 = LayerS(controller.W{1, 3}, controller.b{1, 3}', 'satlin');
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
% step 5: get control output
% step 6: scale control output
% step 7: compute reachable set of the plant with the new control input and
% initial state
% step 8: go back to step 1, ...


lb = [97; 25.2; -3.5];
ub = [97.5; 25.5; -3.0];
init_set = Star(lb, ub); % initial condition of the plant

N = 50; % number of control steps

X0 = init_set; % step 0: initial state of plant
S = X0;


tf_outs = [];
for i=1:N
    % step 1: normalizing state
    norm_X = X0.affineMap(norm_mat, []); % normalized state
    % step 2: compute brake output of rl_controller
    brake = rl_controller.reach(norm_X);
    % step 3: compose brake output and speed
    speed_brakes = get_speed_brakes(brake, norm_X);
    % step 4: compute transformer output
    tf_outs = get_tf_outs(transformer, speed_brakes);
    % step 5: get control output
    controls = get_controls(tf_outs, norm_X);
    % step 6: scale control output
    scaled_controls = scale_controls(controls, scale_mat);
    new_control = Star.get_hypercube_hull(scaled_controls);
    % step 7: compute reachable set of the plant with the new control input and
    % initial state
    X0 = plant.stepReachStar(X0, new_control.toStar);  
    S = [S X0];
end


% plot reachable set
subplot(1,2,1);
Star.plotBoxes_2D_noFill(S, 1, 2, 'b');
xlabel('distance');
ylabel('speed');
title('Speed vs. Distance');
subplot(1,2,2);
Star.plotBoxes_2D_noFill(S, 2, 3, 'b');
xlabel('speed');
ylabel('acceleration');
title('Acceleration vs. Speed');
saveas(gcf, 'reachSet2.pdf');
% save reachable set for computing TTC
save reachSet2.mat S;

function speed_brakes = get_speed_brakes(brake, norm_X)
    n = length(brake);
    speed_brakes = [];
    
    for i=1:n
        V = [norm_X.V(2, :);brake(i).V];
        speed_brakes = [speed_brakes Star(V, brake(i).C, brake(i).d)];
        
    end 
end

function tf_outs = get_tf_outs(transformer, speed_brakes)
    n = length(speed_brakes);
    tf_outs = [];
    parfor i=1:n
        tf_out = transformer.reach(speed_brakes(i));
        tf_outs = [tf_outs tf_out];
    end
end

function controls = get_controls(tf_outs, norm_X)
    n = length(tf_outs);
    controls = [];
    for i=1:n
        V = [norm_X.V(2, :);tf_outs(i).V];
        controls = [controls Star(V, tf_outs(i).C, tf_outs(i).d)];        
    end
end

function scaled_controls = scale_controls(controls, scale_mat)
    n = length(controls);
    scaled_controls = [];
    for i=1:n
        scaled_controls = [scaled_controls controls(i).affineMap(scale_mat, [])];
    end
end