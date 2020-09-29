load controller_7_20.mat;

weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = [];
for i=1:(n-1)
    L = LayerS(weights{1, i}, bias{i, 1}, 'poslin');
    Layers = [Layers L];
end

L = LayerS(weights{1, n}, bias{n, 1}, 'purelin');

Layers = [Layers L];

Controller = FFNNS(Layers); % feedforward neural network controller
reachStep = 0.01; % time step for reachability analysis of the plant
controlPeriod = 0.1; % Ts = 0.1, sampling time for control signal from neural network controller
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Plant = NonLinearODE(6, 1, @car_dynamics, reachStep, controlPeriod, output_mat);

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial condition of the Plant

% x = [x_lead v_lead internal_acc_lead x_ego v_ego internal_acc_ego]'

% initial condition of x_lead
x_lead = cell(10, 1);
x_lead{1, 1} = [108 110];
x_lead{2, 1} = [106 108];
x_lead{3, 1} = [104 106];
x_lead{4, 1} = [102 104];
x_lead{5, 1} = [100 102];
x_lead{6, 1} = [98 100];
x_lead{7, 1} = [96 98];
x_lead{8, 1} = [94 96];
x_lead{9, 1} = [92 94];
x_lead{10, 1} = [90 92];
% initial condition of v_lead
v_lead = [32 32.2];
% initial condition of x_internal_lead
internal_acc_lead = [0 0];
% initial condition of x_ego
x_ego = [10 11]; 
% initial condition of v_ego
v_ego = [30 30.2];
% initial condition of x_internal_ego
internal_acc_ego = [0 0];

n = length(x_lead);

init_set(n) = Star();

for i=1:n
    x1 = x_lead{i, 1};
    lb = [x1(1); v_lead(1); internal_acc_lead(1); x_ego(1); v_ego(1); internal_acc_ego(1)];
    ub = [x1(2); v_lead(2); internal_acc_lead(2); x_ego(2); v_ego(2); internal_acc_ego(2)];
    init_set(i) = Star(lb, ub);
end

% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

lb = [30; 1.4];
ub = [30; 1.4];

input_ref = Star(lb, ub); % input reference set

N = 50;  

n_cores = 4; % number of cores 

% safety specification: relative distance > safe distance
% dis = x_lead - x_ego  
% safe distance between two cars, see here 
% https://www.mathworks.com/help/mpc/examples/design-an-adaptive-cruise-control-system-using-model-predictive-control.html
% dis_safe = D_default + t_gap * v_ego;

t_gap = 1.4;
D_default = 10;

% safety specification: x_lead - x_ego - t_gap * v_ego - D_default > 0
% unsafe region: x_lead - x_ego - t_gap * v_ego <= D_default 

unsafe_mat = [1 0 0 -1 -t_gap 0];
unsafe_vec = [D_default];

safe = zeros(1, n); % safety status
reachTime = zeros(1, n); % reach Time
safetyCheckingTime = zeros(1,n); % safety checking Time
verificationTime = zeros(1, n); % verification tim = reach time + safety checking time
reachPRM.numSteps = N;
reachPRM.reachMethod = 'approx-star';
reachPRM.ref_input = input_ref;
reachPRM.numCores = n_cores;
for i=1:n
    reachPRM.init_set = init_set(i);
    [~, reachTime(i)] = ncs.reach(reachPRM);
    [safe(i), safetyCheckingTime(i)] = ncs.check_safety(unsafe_mat, unsafe_vec, n_cores);
    verificationTime(i) = reachTime(i) + safetyCheckingTime(i);
end

save verify_controller_7_20.mat safe verificationTime;


