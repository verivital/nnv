load controller_5_20.mat;

path_out = [path_results(), filesep, 'ACC', filesep];	

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
% /* car model
Tr = 0.01; % reachability time step for the plant
Tc = 0.1; % control period of the plant
% output matrix
C = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % output matrix
Plant = NonLinearODE(6, 1, @car_dynamics, Tr, Tc, C);

ncs = NonlinearNNCS(Controller, Plant); % the neural network control system

%% Construct initial conditions
% initial condition of the Plant

% x = [x_lead v_lead internal_acc_lead x_ego v_ego internal_acc_ego]'

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
x_ego = [10 11]; 
% initial condition of v_ego
v_ego = [30 30.2];
% initial condition of x_internal_ego
internal_acc_ego = [0 0];

n = length(v_lead);

init_set(n) = Star();

for i=1:n
    lb = [x_lead(1); v_lead{i}(1); internal_acc_lead(1); x_ego(1); v_ego(1); internal_acc_ego(1)];
    ub = [x_lead(2); v_lead{i}(2); internal_acc_lead(2); x_ego(2); v_ego(2); internal_acc_ego(2)];
    init_set(i) = Star(lb, ub);
end

% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

input_ref = [30; 1.4];

%% Verification

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
unsafe_vec = D_default;
unsafeRegion = HalfSpace(unsafe_mat, unsafe_vec);

reachPRM.ref_input = [30; 1.4];
reachPRM.numSteps = 50;
reachPRM.reachMethod = 'approx-star';

if ~exist('numCores')
    reachPRM.numCores = 4;
else
    reachPRM.numCores = numCores;
end

safe = cell(1, n); % safety status
VT = zeros(1,n); % verification time
counterExamples = cell(1,n); % counter examples
for i=1:n
    reachPRM.init_set = init_set(i);
    [safe{i}, counterExamples{i}, VT(i)] = ncs.verify(reachPRM, unsafeRegion);
end

%% Safe verification results
save([path_out, 'verify_nonlinear_ACC.mat'], 'safe', 'VT', 'counterExamples');	


%% Print verification results to screen
fprintf('\n=======================================================');
fprintf('\nVERIFICATION RESULTS FOR ACC WITH NONLINEAR PLANT MODEL');
fprintf('\n=======================================================');
fprintf('\nv_lead(0)                 approx-star    ');
fprintf('\n                       safety  |  VT   ');
for i=1:n
fprintf('\n[%d %d]              %s  |  %3.3f   ', v_lead{i}(1), v_lead{i}(2), safe{i}, VT(i));
end
fprintf('\n-------------------------------------------------------');
fprintf('\nTotal verification time:      %3.3f', sum(VT));

%% Print verification results to a file

fid = fopen([path_out, 'table3_nonlinear_ACC.txt'], 'wt');

fprintf(fid,'\n=======================================================');
fprintf(fid,'\nVERIFICATION RESULTS FOR ACC WITH NONLINEAR PLANT MODEL');
fprintf(fid,'\n=======================================================');
fprintf(fid,'\nv_lead(0)                 approx-star    ');
fprintf(fid,'\n                       safety  |  VT   ');
for i=1:n
fprintf(fid,'\n[%d %d]              %s  |  %3.3f   ', v_lead{i}(1), v_lead{i}(2), safe{i}, VT(i));
end
fprintf(fid,'\n-------------------------------------------------------');
fprintf(fid,'\nTotal verification time:      %3.3f', sum(VT));
fclose(fid);

% %% Plot counter examples	
% cI = counterExamples{1};	
% cI = cell2mat(cI);	
% d_rel = [1 0 0 -1 0 0]*cI;	
% d_safe = [0 0 0 1.4 0 0]*cI + 10;	
% 
% figure; 	
% T = 0:1:50;	
% plot(T, d_rel, 'blue');	
% hold on;	
% plot(T, d_safe, 'red');	
% 
% xlabel('Control Time Steps', 'FontSize', 13);	
% ylabel('Distance', 'FontSize', 13);	
% xticks([0:5:50]);	
% title('Actual Distance (blue) vs. Safe Distance (red)');	
% 
% saveas(gcf, [path_out, 'verify_nonlinear_acc_cex.png']);