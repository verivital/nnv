% Rechability analysis for Emergency Braking System
% [Tran et al., EMSOFT'2019]

%% Load components

% reinforcement learning controller
load('models/controller.mat'); % controller with saturation activation function at the output
rl_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
rl_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
rl_layer3 = LayerS(W{1, 3}, b{1, 3}', 'satlin');
layers_controller = {rl_layer1; rl_layer2; rl_layer3};
rl_controller = NN(layers_controller);

% normalization matrix (plant to controller)
norm_mat = [1/250 0 0; 0 3.6/120 0; 0 0 1/20]; % normalized matrix

% transformation NN
load('models/transform.mat');
tf_layer1 = LayerS(W{1, 1}, b{1, 1}', 'poslin');
tf_layer2 = LayerS(W{1, 2}, b{1, 2}', 'poslin');
tf_layer3 = LayerS(W{1, 3}, b{1, 3}', 'purelin');
layers_transformer = {tf_layer1; tf_layer2; tf_layer3};
transformer = NN(layers_transformer);

% control signal scale (from transformer output to plant input)
scale_mat = [-15*120/3.6 15*120/3.6];

% plant (Discrete linear model)
A = [1 -1/15 0; 0 1 0; 0 0 0];
B = [0;1/15;1];
C = [1 0 0;0 1 0; 0 0 1];
D = [0;0;0];
Ts = 0.1;
plant = DLinearODE(A, B, C, D, Ts);


%% Reachability analysis

% computation steps
% step 0: initial state of plant
% step 1: normalizing state
% step 2: compute brake output of rl_controller
% step 3: compose brake output and normalized speed
% step 4: compute transformer output
% step 5: get control output
% step 6: scale control output
% step 7: compute reachable set of the plant with the new control input and initial state
% step 8: go back to step 1, ...

% Define initial conditions
lb = [97; 25.2; 0];
ub = [97.5; 25.5; 0];
init_set = Star(lb, ub); % initial condition of the plant

N = 90; % number of control steps
% N = 10;

X0 = init_set; % step 0: initial state of plant
S = X0;

start = tic; 
tf_outs = [];
for i=1:N
    % step 1: normalizing state
    norm_X = X0.affineMap(norm_mat, []); % normalized state
    % step 2: the full brake is always 1
    % step 3: compose brake output and speed
    speed_brakes = get_speed_full_brakes(norm_X);
    % step 4: compute transformer output
    tf_outs = get_tf_outs(transformer, speed_brakes);
    % step 5: get control output
    controls = get_controls(tf_outs, norm_X);
    % step 6: scale control output
    scaled_controls = scale_controls(controls, scale_mat);
    % step 7: compute reachable set of the plant with the new control input and initial state
    X0 = stepReachPlant(A, B, X0, scaled_controls); 
    % X0 = plant.stepReachStar(X0, scaled_controls);
    % X0 = Star.get_hypercube_hull(X0);
    % X0 = X0.toStar;
    S = [S X0];
    % step 8: start again for next time step...
end

reachTime = toc(start);

%% Visualization

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
%saveas(gcf, 'reachSet_fullbrake.pdf');
%save reachable set for computing TTC
% save reachSet.mat S;



%% Helper Functions

% Get controller output (always 1 for full brake)
function speed_brakes = get_speed_full_brakes(norm_X) 
    V1 = zeros(1, norm_X.nVar + 1);
    V1(1,1) = 1;
    new_V = [norm_X.V(2,:); V1];
    speed_brakes = Star(new_V, norm_X.C, norm_X.d);
    speed_brakes.predicate_lb = -1*ones(speed_brakes.nVar,1);
    speed_brakes.predicate_ub = ones(speed_brakes.nVar,1);
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
        controls = [controls Star(V, tf_outs(i).C, tf_outs(i).d)];        
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