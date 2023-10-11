% This script performs rechability analysis for Emergency Braking System
% author: Dung Tran
% date: 6/3/2019

%% Show AEBS scenario
imshow("AEBS_scenario.png");
% From Tran et al., "Safety Verification of Cyber-Physical Systems with
% Reinforcement Learning Control", EMSOFT'19

%% Begin verification
% normalization 
norm_mat = [1/250 0 0; 0 3.6/120 0; 0 0 1/20]; % normalized matrix

% reinforcement learning controller
load('models/controller.mat'); % controller with saturation activation function at the output
rl_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
rl_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
rl_layer3 = LayerS(W{1, 3}, b{1, 3}', 'satlin');
layers_controller = {rl_layer1; rl_layer2; rl_layer3};
rl_controller = NN(layers_controller);

% transformation 
load('models/transform.mat');
tf_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
tf_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
tf_layer3 = LayerS(W{1, 3}, b{1, 3}', 'purelin');
layers_transformer = {tf_layer1; tf_layer2; tf_layer3};
transformer = NN(layers_transformer);

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
init_set = Star(lb, ub); % initial condition of the plant

N = 50; % number of control steps

X0 = init_set; % step 0: initial state of plant
S = X0;

start = tic; 
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
    % step 7: compute reachable set of the plant with the new control input and
    % initial state
    X0 = stepReachPlant(A, B, X0, scaled_controls);  
    S = [S X0];
end
% save reachSet.mat S;

reachTime = toc(start);

% plot reachable set
subplot(1,2,1);
Star.plotBoxes_2D_noFill(S, 1, 2, 'b');
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
set(gca,'FontSize',16)
%title('Speed vs. Distance');
subplot(1,2,2);
Star.plotBoxes_2D_noFill(S, 2, 3, 'b');
xlabel('Velocity (m/s)');
ylabel('Acceleration (m/s^2)');
set(gca,'FontSize',16)
% saveas(gcf, 'reachSet_approx.pdf');
%save reachable set for computing TTC
% save reachSet.mat S;


%% Helper Functions

% Get the inputs to the transformer
function speed_brakes = get_speed_brakes(brake, norm_X)
    n = length(brake);
    speed_brakes = [];
    for i=1:n
        V1 = zeros(size(brake(i).V));
        V1(1:length(norm_X.V(2, :))) = norm_X.V(2, :);
        V = [V1; brake(i).V];
%         V = [norm_X.V(2, :); brake(i).V];
        speed_brakes = [speed_brakes Star(V, brake(i).C, brake(i).d, brake(i).predicate_lb, brake(i).predicate_ub)];
    end 
end

% Compute transformer reachability
function tf_outs = get_tf_outs(transformer, speed_brakes)
    n = length(speed_brakes);
    tf_outs = [];
    for i=1:n
        tf_out = transformer.reach(speed_brakes(i)); % approx-star (default)
        tf_outs = [tf_outs tf_out];
    end
end

% Get control values
function controls = get_controls(tf_outs, norm_X)
    n = length(tf_outs);
    controls = [];
    for i=1:n
        V1 = zeros(size(tf_outs(i).V));
        V1(1:length(norm_X.V(2, :))) = norm_X.V(2, :);
        V = [V1; tf_outs(i).V];
        controls = [controls Star(V, tf_outs(i).C, tf_outs(i).d, tf_outs(i).predicate_lb, tf_outs(i).predicate_ub)];        
    end
    % in1 = norm_X.affineMap([0 1 0],[]); % get speed set
    % for i=1:n
    %     cSet = in1.concatenate(tf_outs(i));
    %     controls = [controls cSet];        
    % end
end

% Normalized control output
function scaled_controls = scale_controls(controls, scale_mat)
    n = length(controls);
    scaled_controls = [];
    for i=1:n
        scaled_controls = [scaled_controls controls(i).affineMap(scale_mat, [])];
    end
end

function X = stepReachPlant(A, B, X0, scaled_controls)
    n = length(scaled_controls);
    X = X0.affineMap(A, []);
    U = [];
    for i=1:n
        U = [U scaled_controls(i).affineMap(B, [])];
    end
    next_X = [];
    for i=1:n
        V1 = zeros(size(U(i).V));
        V1(1:size(X.V,1), 1:size(X.V,2)) = X.V;
        V = V1 + U(i).V;
        next_X = [next_X Star(V, U(i).C, U(i).d, -1*ones(U(i).nVar,1), ones(U(i).nVar,1))];
    end
    X = Star.get_hypercube_hull(next_X);
    X = X.toStar;
end