% Reachability analysis for Discrete Linear ACC model
% Dung Tran: 9/30/2019



%% System model
% x1 = lead_car position
% x2 = lead_car velocity
% x3 = lead_car internal state

% x4 = ego_car position
% x5 = ego_car velocity
% x6 = ego_car internal state

% lead_car dynamics
%dx(1,1)=x(2);
%dx(2,1) = x(3);
%dx(3,1) = -2 * x(3) + 2 * a_lead - mu*x(2)^2;

% ego car dynamics
%dx(4,1)= x(5); 
%dx(5,1) = x(6);
%dx(6,1) = -2 * x(6) + 2 * a_ego - mu*x(5)^2;


% let x7 = -2*x(3) + 2 * a_lead -> x7(0) = -2*x(3)(0) + 2*alead
% -> dx7 = -2dx3


A = [0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 0 0 0 1; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 -2 0; 0 0 -2 0 0 0 0];
B = [0; 0; 0; 0; 0; 2; 0];
C = [1 0 0 -1 0 0 0; 0 1 0 0 -1 0 0; 0 0 0 0 1 0 0];  % feedback relative distance, relative velocity, longitudinal velocity
D = [0; 0; 0]; 

plant = LinearODE(A, B, C, D); % continuous plant model

plantd = plant.c2d(0.1); % discrete plant model

% the neural network provides a_ego control input to the plant
% a_lead = -5 


%% Controller
load controller_5_20.mat;

weights = network.weights;
bias = network.bias;
n = length(weights);
n = length(weights);
Layers = [];
for i=1:n - 1
    L = LayerS(weights{1, i}, bias{i, 1}, 'poslin');
    Layers = [Layers L];
end
L = LayerS(weights{1, n}, bias{n, 1}, 'purelin');
Layers = [Layers L];
Controller = FFNNS(Layers); % feedforward neural network controller


%% NNCS 

ncs = DLinearNNCS(Controller, plantd); % a discrete linear neural network control system


%% Initial Set of states and reference inputs

% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

ref_input = [30; 1.4];

% initial condition of the plant

% initial position of lead car x_lead
x_lead = [90 92];
% initial condition of v_lead
v_lead = cell(10, 1);
v_lead{1, 1} = [29 30];
v_lead{2, 1} = [28 29];
v_lead{3, 1} = [27 28];
v_lead{4, 1} = [26 27];
v_lead{5, 1} = [25 26];
v_lead{6, 1} = [24 25];
v_lead{7, 1} = [23 24];
v_lead{8, 1} = [22 23];
v_lead{9, 1} = [21 22];
v_lead{10, 1} = [20 21];

% initial condition of x_internal_lead
internal_acc_lead = [0 0];
% initial condition of x_ego
x_ego = [30 31]; 
% initial condition of v_ego
v_ego = [30 30.5];
% initial condition of x_internal_ego
internal_acc_ego = [0 0];

a_lead = -5;
% initial condition for new introduced variable 
x7_0 = [2*a_lead 2*a_lead]; % x7 = -2*x(3) + 2 * a_lead -> x7(0) = -2*x(3)(0) + 2*alead = -2*0 + 2*-2 = -4

N = length(v_lead);
for i=1:N
    lb = [x_lead(1); v_lead{i}(1); internal_acc_lead(1); x_ego(1); v_ego(1); internal_acc_ego(1); x7_0(1)];
    ub = [x_lead(2); v_lead{i}(2); internal_acc_lead(2); x_ego(2); v_ego(2); internal_acc_ego(2); x7_0(2)];
    init_set(i) = Star(lb, ub);
end



%% Reachability Analysis && Verification
numSteps = 50; 
numCores = 4; 
% safety property: actual distance > alpha * safe distance <=> d = x1 - x4  > alpha * d_safe = alpha * (1.4 * v_ego + 10)

% usafe region: x1 - x4 <= alpha * (1.4 * v_ego + 10)
alpha = 1;
unsafe_mat = [1 0 0 -1 -alpha*1.4 0 0];
unsafe_vec = alpha*10; 

% N = 1; % just for testing
safe_exact = cell(1,N);
VT_exact = zeros(1, N);
counterExamples_exact = cell(1, N);
dis_ACC_exact = cell(1,N);

for i=1:N
    [safe_exact{i}, counterExamples_exact{i}, VT_exact(i)] = ncs.verify(init_set(i), ref_input, numSteps, 'exact-star', numCores, unsafe_mat, unsafe_vec); 
    dis_ACC_exact{i} = ncs; %store for plotting reachable sets
end


safe_approx = cell(1, N);
VT_approx = zeros(1, N);
counterExamples_approx = cell(1, N);
dis_ACC_approx = cell(1, N);

for i=1:N
    [safe_approx{i}, counterExamples_approx{i}, VT_approx(i)] = ncs.verify(init_set(i), ref_input, numSteps, 'approx-star', numCores, unsafe_mat, unsafe_vec);
    dis_ACC_approx{i} = ncs; % store for plotting reachable sets
end


%% Safe verification results

save discrete_ACC_verification_results.mat safe_exact VT_exact counterExamples_exact safe_approx VT_approx counterExamples_approx dis_ACC_exact dis_ACC_approx;

%% Print verification results to screen
fprintf('\n======================================================');
fprintf('\nVERIFICATION RESULTS FOR ACC WITH DISCRETE PLANT MODEL');
fprintf('\n======================================================');
fprintf('\nv_lead(0)        exact-star          approx-star    ');
fprintf('\n               safety  |  VT        safety  |  VT   ');
fprintf('\n[%d %d]       %s  |  %3.3f        %s  |  %3.3f   ', v_lead{1}(1), v_lead{1}(2), safe_exact{1}, VT_exact(1), safe_approx{1}, VT_approx(1));
for i=2:N
fprintf('\n[%d %d]        %s  |  %3.3f        %s  |  %3.3f   ', v_lead{i}(1), v_lead{i}(2), safe_exact{i}, VT_exact(i), safe_approx{i}, VT_approx(i));
end

%% Print verification results to a file

fid = fopen('discrete_ACC_verification_results.txt', 'wt');
fprintf(fid,'\n======================================================');
fprintf(fid,'\nVERIFICATION RESULTS FOR ACC WITH DISCRETE PLANT MODEL');
fprintf(fid,'\n======================================================');
fprintf(fid,'\nx_lead(0)        exact-star          approx-star    ');
fprintf(fid,'\n               safety  |  VT        safety  |  VT   ');
fprintf(fid,'\n[%d %d]       %s  |  %3.3f        %s  |  %3.3f   ', v_lead{1}(1), v_lead{1}(2), safe_exact{1}, VT_exact(1), safe_approx{1}, VT_approx(1));
for i=2:N
fprintf(fid,'\n[%d %d]        %s  |  %3.3f        %s  |  %3.3f   ', v_lead{i}(1), v_lead{i}(2), safe_exact{i}, VT_exact(i), safe_approx{i}, VT_approx(i));
end
fclose(fid);

%% END
