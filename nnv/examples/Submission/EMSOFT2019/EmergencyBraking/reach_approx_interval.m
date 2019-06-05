% This script performs rechability analysis for Emergency Braking System
% author: Dung Tran
% date: 6/3/2019

% normalization 
norm_mat = [1/250 0 0; 0 3.6/120 0; 0 0 1/20]; % normalized matrix

% reinforcement learning controller
load controller.mat; % controller with saturation activation function at the output
rl_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
rl_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
rl_layer3 = LayerS(W{1, 3}, b{1, 3}', 'satlin');
rl_controller = FFNNS([rl_layer1 rl_layer2 rl_layer3]); 

% transformation 
load transform.mat;
tf_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
tf_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
tf_layer3 = LayerS(W{1, 3}, b{1, 3}', 'purelin');
transformer = FFNNS([tf_layer1 tf_layer2 tf_layer3]);

% control signal scale
scale_mat = [-15*120/3.6 15*120/3.6];

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


lb = [97; 25.2; 0];
ub = [97.5; 25.5; 0];


%vmax = 4;
%lb = [10; vmax; 0]; % using for searching
%ub = [20; vmax; 0];
%init_set = Star(lb, ub); % initial condition of the plant
init_set = Polyhedron('lb', lb, 'ub', ub);
N = 2; % number of control steps

X0 = init_set; % step 0: initial state of plant
S = cell(1, N+1);
S{1,1} = X0;

start = tic; 
tf_outs = [];
for i=1:N
    % step 1: normalizing state
    norm_X = X0.affineMap(norm_mat); % normalized state
    % step 2: compute brake output of rl_controller
    brake = rl_controller.reach(norm_X, 'exact-polyhedron',4);
    % step 3: compose brake output and speed
    speed_brakes = get_speed_brakes(brake, norm_X);
    % step 4: compute transformer output
    tf_outs = get_tf_outs(transformer, speed_brakes);
    % step 5: get control output
    controls = get_controls(tf_outs, norm_X);
    % step 6: scale control output
    scaled_controls = scale_controls(controls, scale_mat);
    % step 7: compute reachable set of the plant with the new control input and
    % initial state
    [X0, Y] = stepReachPlant(A, B, X0, scaled_controls);  
    S{1, i+1} = Y;
end
save reachSet_Interval.mat S;

reachTime = toc(start);
% plot reachable set
% plot reachable set

figure;
for i=1:N+1
    Si = S{1, i}.affineMap([1 0 0; 0 1 0]);
    Si.plot;
    hold on;
end
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
set(gca,'FontSize',16)
title('Speed vs. Distance');

function speed_brakes = get_speed_brakes(brake, norm_X)
    
    B1 = Reduction.hypercubeHull(brake);
    norm_X.outerApprox;
    speed = norm_X.affineMap([0 1 0]);
    
    speed_brakes = [speed B1.toPolyhedron];  
    speed_brakes = Conversion.concatenatePolyhedron(speed_brakes);
    
end

function tf_outs = get_tf_outs(transformer, speed_brakes)
    n = length(speed_brakes);
    tf_outs = [];
    parfor i=1:n
        tf_out = transformer.reach(speed_brakes(i), 'exact-polyhedron');
        tf_outs = [tf_outs tf_out];
    end
end

function controls = get_controls(tf_outs, norm_X)
    n = length(tf_outs);
    controls = [];
    speed = norm_X.affineMap([0 1 0]);
    parfor i=1:n
        Ps = [speed tf_outs(i)];
        controls = [controls Conversion.concatenatePolyhedron(Ps)];        
    end
end

function scaled_controls = scale_controls(controls, scale_mat)
    n = length(controls);
    scaled_controls = [];
    parfor i=1:n
        scaled_controls = [scaled_controls controls(i).affineMap(scale_mat)];
    end
    scaled_controls = Reduction.hypercubeHull(scaled_controls);
    scaled_controls = scaled_controls.toPolyhedron;
end

function [X, Y] = stepReachPlant(A, B, X0, scaled_controls)
    
    n = length(scaled_controls);
    X = X0.affineMap(A);
    U = [];
    parfor i=1:n
        U = [U scaled_controls(i).affineMap(B)];
    end
    
    next_X = [];
    parfor i=1:n
        V = X + U(i);
        next_X = [next_X V];
    end
    
    
    X = Reduction.hypercubeHull(next_X);
    X = X.toPolyhedron;
    Y = X;
end