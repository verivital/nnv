load controller_3_20.mat;

n = length(weights);
Layers = [];
for i=1:n - 1
    L = Layer(weights{1, i}, bias{i, 1}, 'ReLU');
    Layers = [Layers L];
end

L = Layer(weights{1, n}, bias{n, 1}, 'Linear');

Layers = [Layers L];

Controller = FFNN(Layers); % feedforward neural network controller

pwd

[ha_cfg, ha_name] = createConfigFromSpaceEx('acc.xml', 'acc.cfg');

ha_name
ha_cfg

% TODO: for now, assuming there is 1 network component and 1 child
% automaton, generalize later, just want to pull out the ODEs for now
ha = ha_cfg.root.template.children.entrySet().iterator().next().getValue().child;

% TODO: also assume 1 mode, generalize later
mode = ha.modes.entrySet().iterator().next().getValue();

% TODO: iterator over all variables, get flows for each
flows = mode.flowDynamics.entrySet().iterator();

vnames = ha.variables.iterator();

while vnames.hasNext()
    syms(vnames.next());
end

constants = ha.constants.entrySet().iterator();

while constants.hasNext()
    syms(constants.next().getKey());
end

syms x; % todo: pull in vectorized prefix name

f_ode_str = '';

iflow = 1;
while flows.hasNext()
    f_next = flows.next();
    vname = f_next.getKey()
    f = char(f_next.getValue().toDefaultString())
    
    % TODO: generalize or add a renaming variable pass, this assumes all
    % variables start with names like xN where N is a possibly multi-digit
    % number
    f_matlab = regexprep(f, 'x[0-9]*', 'x(${strrep($0, ''x'', '''')})')
    
    f_ode = strcat(vname, ' = ', f);
    
    %f_ode_sym(iflow) = str2sym(f_ode)
    %f_ode_sym{iflow} = f;
    f_ode_sym(iflow) = str2sym(f_matlab)
    fh_sym(iflow) = symfun(f_ode_sym(iflow), [x1 x2 x3 x4 x5 x6 a_ego a_lead mu])
    
    funcstr{iflow} = str2func(char(f_ode_sym(iflow)));
    
    % char(13) = '\r'
    f_ode_str = strcat(f_ode_str, newline, char(13), 'dx(', num2str(iflow), ', 1) = ', f_matlab, ';');

%function [dx]=acc_St1_FlowEq(t,x,u)
%dx(1,1) = x(2);
%dx(2,1) = x(3);
%dx(3,1) = 2*u(1) - 2*x(3) - u(3)*x(2)^2;
%dx(4,1) = x(5);
%dx(5,1) = x(6);
%dx(6,1) = 2*u(2) - 2*x(6) - u(3)*x(5)^2;
    
    
    iflow = iflow + 1;
end

%dx(1,1)=x(2);
%dx(2,1) = x(3);
%dx(3,1) = -2 * x(3) + 2 * a_lead - mu*x(2)^2;
% ego car dyanmics
%dx(4,1)= x(5); 
%dx(5,1) = x(6);
%dx(6,1) = -2 * x(6) + 2 * a_ego - mu*x(5)^2;

%Plant = NonLinearODE(6, 1, @car_dynamics);

%Plant = NonLinearODE(6, 1, symfun(f_ode_sym, [x1 x2 x3 x4 x5 x6 a_ego a_lead mu]));

printDynamicsFile(pwd,'acc_dyn',f_ode_str);
copyfile( strcat(pwd,filesep,'auxiliary', filesep, 'acc_dyn.m'), strcat(pwd, filesep));

Plant = NonLinearODE(6, 1, @acc_dyn);

%return

Plant.set_timeStep(0.01); % time step for reachability analysis of the plant
Plant.set_tFinal(0.1); % Ts = 0.1, sampling time for control signal from neural network controller
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

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
for i=1:n
    [~, reachTime(i)] = ncs.reach('approx-star', init_set(i), input_ref, n_cores, N);
    [safe(i), safetyCheckingTime(i)] = ncs.check_safety(unsafe_mat, unsafe_vec, n_cores);
    verificationTime(i) = reachTime(i) + safetyCheckingTime(i);
end

save verify_controller_3_20.mat safe verificationTime;


