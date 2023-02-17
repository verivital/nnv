% Reachability analysis for Discrete Linear ACC model
% Dung Tran: 9/30/2019

path_out = [path_results(), filesep, 'ACC', filesep];	

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
v_lead = cell(6, 1);
v_lead{1, 1} = [29 30];
v_lead{2, 1} = [28 29];
v_lead{3, 1} = [27 28];
v_lead{4, 1} = [26 27];
v_lead{5, 1} = [25 26];
v_lead{6, 1} = [24 25];

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
numSteps =10;
% safety property: actual distance > alpha * safe distance <=> d = x1 - x4  > alpha * d_safe = alpha * (1.4 * v_ego + 10)

% usafe region: x1 - x4 <= alpha * (1.4 * v_ego + 10)
alpha = 1;
unsafe_mat = [1 0 0 -1 -alpha*1.4 0 0];
unsafe_vec = alpha*10;

% reachability parameters

reachPRM.ref_input = ref_input;
reachPRM.numSteps = numSteps;
reachPRM.reachMethod = 'approx-star';

if ~exist('numCores')
    reachPRM.numCores = 4;
else
    reachPRM.numCores = numCores;
end


unsafeRegion = HalfSpace(unsafe_mat, unsafe_vec);

% N = 1; % just for testing

safe_approx = cell(1, N);
VT_approx = zeros(1, N);
counterExamples_approx = cell(1, N);
dis_ACC_approx = cell(1, N);

for i=1:N
    reachPRM.init_set = init_set(i);
    [safe_approx{i}, counterExamples_approx{i}, VT_approx(i)] = ncs.verify(reachPRM, unsafeRegion);
end


%% Safe verification results
save([path_out, 'linear_ACC.mat'], 'safe_approx', 'VT_approx', 'counterExamples_approx');

%% Print verification results to screen
fprintf('\n======================================================');
fprintf('\nVERIFICATION RESULTS FOR ACC WITH DISCRETE PLANT MODEL');
fprintf('\n======================================================');
fprintf('\nv_lead(0)                 approx-star    ');
fprintf('\n                       safety  |  VT   ');
for i=1:N
fprintf('\n[%d %d]                %s    %3.3f   ', v_lead{i}(1), v_lead{i}(2), safe_approx{i}, VT_approx(i));
end
fprintf('\n-------------------------------------------------------');
fprintf('\nTotal verification time:      %3.3f', sum(VT_approx));

%% Print verification results to a file

fid = fopen([path_out, 'table3_linear_ACC.txt'], 'wt');	

fprintf(fid,'\n======================================================');
fprintf(fid,'\nVERIFICATION RESULTS FOR ACC WITH DISCRETE PLANT MODEL');
fprintf(fid,'\n======================================================');
fprintf(fid,'\nv_lead(0)                 approx-star    ');
fprintf(fid,'\n                       safety  |  VT   ');
for i=1:N
fprintf(fid,'\n[%d %d]                %s  |  %3.3f   ', v_lead{i}(1), v_lead{i}(2), safe_approx{i}, VT_approx(i));
end
fprintf(fid,'\n-------------------------------------------------------');
fprintf(fid,'\nTotal verification time:      %3.3f', sum(VT_approx));
fclose(fid);

%% END
