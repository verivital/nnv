%% NNCS reachability exercise of a discrete linearized version of the ACC

% Load controller
load('controller.mat');

% Load plant dynamics (discrete linear)
load("sysd.mat");
Plant = DLinearODE(sysd.A, sysd.B, sysd.C, sysd.D, sysd.Ts);

% Create neural network control system
ncs = NNCS(); % create NNCS with Controller + Plant

% initial condition of the system

% state variables
% x = [x_lead; v_lead; x_ego; v_ego]'

% Requirements for initial set
% Positions (x_lead, x_ego)
%   - Initial distance between cars must be 60 and 100
%   - Size must be <= 15, i.e., (upper_bounds - lower_bound of x) <= 15
% Velocity (v_lead, v_ego)
%   - Initial speeds must be between 30 and 35
%   - Size must be <= 0.5

% define bounds as vectors
lb = ;
ub = ;

init_set = Star(lb, ub); % define init_set a Star with upper and lower bounds

% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

lb = [1.4; 30];
ub = [1.4; 31];
input_ref = Star(lb, ub);

% define reachability parameters for NNCS
reachPRM.ref_input = input_ref;
reachPRM.numCores = 1;
reachPRM.numSteps = ; % choose number of control steps (please choose a value less than 50) 
reachPRM.reachMethod = 'approx-star';
reachPRM.; % add the initial set to the reachability parameters

[R, reachTime] = ncs.reach(reachPRM); %complete function to do reachability of an NNCS

% plot output relative distance, ego_car velocity and safe_distance
% distance between two cars
% dis = x_lead - x_ego 

% safe distance between two cars, see here 
% https://www.mathworks.com/help/mpc/examples/design-an-adaptive-cruise-control-system-using-model-predictive-control.html

% dis_safe = D_default + t_gap * v_ego;
% the safety specification is: dis >= 0.9 * dis_safe

t_gap = 1.4;
D_default = 10;

% safety specification: distance > alp * safe_distance

map = [1 0 -1 0]; % get distance between two cars and velocity of ego car

N = length(R); % number of output sets

S(N) = Star;
for i=1:N
    S(i) = R(i).affineMap(map, []);
end


% plot safe_distance vs. velocity
alp = 1;
map1 = [0 0 0 alp*t_gap]; % safe distance and velocity of ego car
S1(N) = Star;
for i=1:N
    S1(i) = R(i).affineMap(map1, alp*D_default);
end

controlSteps = 0:1:N-1; % get the control steps for plotting

% Plot sets over time (control steps)
figure;
% Actual distance between cars
Star.plotRanges_2D(S, 1, controlSteps,'b'); % plot sets S, dimensions 1 
ylabel('Safe distances (red) and actual distances (blue)');
xlabel('Control steps');
% Safety distance
hold on;
Star.plotRanges_2D(S1, 1, controlSteps, 'r'); % plot sets S1, dimensions 1 


% Visualize the speeds of the car along control Steps

% Step 1 - Project the Star set (R) into a 1D set for velocity of lead car
Vlead(N) = Star; % initialize set
for i=1:N
    Vlead(i) = R.; % finish function to project 4D set into 1D
end

% Step 2 - Project the Star set (R) into a 1D set for velocity of ego car
Vego(N) = Star; % initialize set
for i=1:N
    Vego(i) = R.; % finish function to project 4D set into 1D
end

% Step 3 - Plot velocity sets over time (control steps)
figure;
% lead car velocity
Star.plotRanges_2D(); % finish function to plot Vlead over the control steps
ylabel('ego car velocity (red) and lead car velocity (blue)');
xlabel('Control steps');
% Safety distance
hold on;
Star.plotRanges_2D(); % % finish function to plot Vego over the control steps
