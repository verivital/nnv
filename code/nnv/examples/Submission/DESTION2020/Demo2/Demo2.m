%% Using NNV to verify a NNCS (ACC)

%% Step 1. Load components

% Plants

% Linear model
A = [0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 0 0 0 1; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 -2 0; 0 0 -2 0 0 0 0];
B = [0; 0; 0; 0; 0; 2; 0];
C = [1 0 0 -1 0 0 0; 0 1 0 0 -1 0 0; 0 0 0 0 1 0 0];  % feedback relative distance, relative velocity, longitudinal velocity
D = [0; 0; 0]; 
acc_linear = LinearODE(A, B, C, D); 

% Nonlinear model
Tr = 0.01; % reachability time step for the plant
Tc = 0.1; % control period of the plant
OutMat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % output matrix
acc_nonlinear = NonLinearODE(6, 1, @car_dynamics, Tr, Tc, OutMat); % nonlinear

% Controller
NNcont = Load_nn('ACCcontroller5_20.mat');

% NNCS 
ncs_nonlinear = NonlinearNNCS(NNcont,acc_nonlinear);
ncs_linear = LinearNNCS(NNcont,acc_linear);

%% Step 2. Set initial conditions (linear)

% Initial states
x_lead = [90 92];       % x1
v_lead = [29 30];       % x2
a_lead = [0 0];           % x3
x_ego = [30 31];        % x4
v_ego = [30 30.2];     % x5
a_ego = [0 0];             % x6
v_new = [-10 -10];     % x7
% Create star set
lb = [x_lead(1);v_lead(1);a_lead(1);x_ego(1);v_ego(1);a_ego(1);v_new(1)];
ub = [x_lead(2);v_lead(2);a_lead(2);x_ego(2);v_ego(2);a_ego(2); v_new(2)];
init_set = Star(lb,ub);
% Reference input for controller
input_ref = [30; 1.4]; % t_gap; v_set

%% Live Reachability Analysis (linear)
numSteps = 50; 
method = 'approx-star';
numCores = 4;
plantNumOfSimSteps = 20;
% plot on-the-fly the distance between two cars x1-x4 and the safe distance d_safe = D_default + t_gap * v_ego = 10 + 1.4 * x(5)
output_mat = [1 0 0 -1 0 0 0]; % plot on-the-fly the distance between two cars x1 - x4 
output_vec = [];
boundary_mat = [0 0 0 0 1.4 0 0]; % plot on-the-fly the safe distance 
boundary_vec = 10;
figureTitle = 'Actual distance (blue) vs. safe distance (red)';

% safe scenarios
figure;
ncs_linear.reachLive(init_set, input_ref, numSteps, 'outputMatrix', output_mat, 'outputVector', output_vec, 'outputSetColor', 'blue', 'boundaryMatrix', boundary_mat, 'boundaryVector', boundary_vec, 'plantNumOfSimSteps', plantNumOfSimSteps, 'figureTitle', figureTitle, 'videoName', 'reachLive_video', 'figureXLabel', 'Time (seconds)', 'figureYLabel', 'Distance (m)', 'videoFrameRate', 40);

%% Step 2. Set initial conditions (nonlinear)

% Initial states
x_lead = [90 92];       % x1
v_lead = [29 30];       % x2
a_lead = [0 0];           % x3
x_ego = [30 31];        % x4
v_ego = [30 30.2];     % x5
a_ego = [0 0];             % x6
% Create star set
lb = [x_lead(1);v_lead(1);a_lead(1);x_ego(1);v_ego(1);a_ego(1)];
ub = [x_lead(2);v_lead(2);a_lead(2);x_ego(2);v_ego(2);a_ego(2)];
init_set = Star(lb,ub);
% Reference input for controller
input_ref = [30; 1.4]; % t_gap; v_set

%% Reachability anaysis (nonlinear)
t_gap = 1.4;
D_default = 10;

% safety specification: x_lead - x_ego - t_gap * v_ego - D_default > 0
% unsafe region: x_lead - x_ego - t_gap * v_ego <= D_default 

unsafe_mat = [1 0 0 -1 -t_gap 0];
unsafe_vec = D_default;
unsafeRegion = HalfSpace(unsafe_mat, unsafe_vec);
% Reachability options
reachPRM.ref_input = [30; 1.4];
reachPRM.numSteps = 50;
reachPRM.reachMethod = 'approx-star';
reachPRM.numCores = 4;
reachPRM.init_set = init_set;

[safe, counterExs, VT] = ncs_nonlinear.verify(reachPRM, unsafeRegion);

%% Print verification results to screen (one case)
fprintf('\n=======================================================');
fprintf('\nVERIFICATION RESULTS FOR ACC WITH NONLINEAR PLANT MODEL');
fprintf('\n=======================================================');
fprintf('\nv_lead(0)                 approx-star    ');
fprintf('\n                                safety  |  VT   ');
fprintf('\n[%d %d]              %s  |  %3.3f   ', v_lead(1), v_lead(2), safe, VT);


%% Print verification results to screen (multiple initial sets)
load('nonlinearACC_results.mat');
n = length(v_lead);
% Difference ->  x_ego = [10 11];      % x4
fprintf('\n=======================================================');
fprintf('\nVERIFICATION RESULTS FOR ACC WITH NONLINEAR PLANT MODEL');
fprintf('\n=======================================================');
fprintf('\nv_lead(0)                approx-star    ');
fprintf('\n                              safety  |  VT   ');
for i=1:n
    fprintf('\n[%d %d]              %s  |  %3.3f   ', v_lead{i}(1), v_lead{i}(2), safe{i}, VT(i));
end
fprintf('\n-------------------------------------------------------');
fprintf('\nTotal verification time:      %3.3f', sum(VT));

%% Plot counter example
cI = counterExs{1,1};	
cI = cell2mat(cI);	
d_rel = [1 0 0 -1 0 0]*cI;	
d_safe = [0 0 0 1.4 0 0]*cI + 10;	

figure; 	
T = 0:1:50;	
plot(T, d_rel, 'blue');	
hold on;	
plot(T, d_safe, 'red');	

xlabel('Control Time Steps', 'FontSize', 13);	
ylabel('Distance', 'FontSize', 13);	
xticks([0:5:50]);	
title('Actual Distance (blue) vs. Safe Distance (red)');	

saveas(gcf, 'counterEx_nonlinear.png');

