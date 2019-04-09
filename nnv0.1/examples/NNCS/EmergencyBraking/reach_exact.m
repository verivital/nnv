% This script performs rechability analysis for Emergency Braking System
% with the RL controller that has saturation activation function
% author: Dung Tran
% date: 6/3/2019

% normalization 
norm_mat = [1/250 0 0; 0 3.6/120 0; 0 0 1/20]; % normalized matrix

% reinforcement learning controller
%load ControllerSat.mat;
load controller.mat; 
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
init_set = Star(lb, ub); % initial condition of the plant

N = 50; % number of control steps, takes arround 800 seconds to compute

X_cell = cell(1, N + 1);
X_cell{1, 1} = init_set;

start = tic;

for i=1:N
    X0 = X_cell{1, i};
    M = length(X0);
    X = [];
    % exact reach for all initial sets X0
    parfor j=1:M        
        % exact reach for a single X0
        X = [X reach_exact_single_X0(rl_controller, transformer, A, B, X0(j), norm_mat, scale_mat)];
    end
    
    X_cell{1, i+1} = X;
end

reachTime = toc(start);

start = tic;
N = length(X_cell);
S = [];
for i=1:N
    fprintf('\nget hypercube hull of the %d^th step reach sets', i);
    B = Star.get_hypercube_hull(X_cell{1, i});
    if isempty(B)
        fprintf('\ni = %d -> empty hull', i);
    else
        S = [S Star.get_hypercube_hull(X_cell{1, i})];
    end
end
save exactReachSet.mat S;
get_ranges_time = toc(start); 

% plot reachable set
figure;
subplot(1,2,1);
Box.plotBoxes_2D_noFill(S, 1, 2, 'b');
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
set(gca,'FontSize',16)
%title('Speed vs. Distance');
subplot(1,2,2);
Box.plotBoxes_2D_noFill(S, 2, 3, 'b');
xlabel('Velocity (m/s)');
ylabel('Acceleration (m/s^2)');
set(gca,'FontSize',16);
%title('Acceleration vs. Speed');
% save reachable set for computing TTC


% plot number of stars of the reachable set over time
N = length(X_cell);
numOfStar = zeros(N, 1);
for i=1:N
    numOfStar(i) = length(X_cell{1, i});
end
figure;
T = 1:1:N;
plot(T, numOfStar, '-x');
xlabel('Time steps');
ylabel('Number of stars');
set(gca,'FontSize',16)
saveas(gcf, 'numOfStar_exactMethod.pdf');


function next_X = reach_exact_single_X0(rl_controller, transformer, A, B, X0, norm_mat, scale_mat)
    
    norm_X = X0.affineMap(norm_mat, []); % normalize state
    brake = rl_controller.reach(norm_X); % compute exact brake
    m = length(brake);
    speed_brakes = []; % get exact speed_brakes
    for j=1:m
         V = [norm_X.V(2, :);brake(j).V];
         speed_brake = Star(V, brake(j).C, brake(j).d);
         speed_brakes = [speed_brakes speed_brake];
    end

    tf_outs = get_tf_outs(transformer, speed_brakes); % compute exact tf_outs
    control = get_controls(tf_outs, norm_X); % get exact control
    scaled_control = scale_controls(control, scale_mat); % scale exact control
    next_X= stepReachPlant(A, B, X0, scaled_control); % compute exact state for next step
        
end

function tf_outs = get_tf_outs(transformer, speed_brakes)
    n = length(speed_brakes);
    tf_outs = [];
    for i=1:n
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

function next_X = stepReachPlant(A, B, X0, scaled_controls)

    n = length(scaled_controls);
    X = X0.affineMap(A, []);
    U = [];
    for i=1:n
        U = [U scaled_controls(i).affineMap(B, [])];
    end
    
    next_X = [];
    for i=1:n
        V = X.V + U(i).V;
        next_X = [next_X Star(V, U(i).C, U(i).d)];
    end
    
end