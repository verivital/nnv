% example call:
% close all ; clear all ; pReach = figure ; pRuntime = figure ; reach_nncACCsystem_video

load controller.mat;
load sysd.mat;

n = length(weights);
Layers = [];
for i=1:n - 1
    L = Layer(weights{1, i}, bias{i, 1}, 'ReLU');
    Layers = [Layers L];
end

L = Layer(weights{1, n}, bias{n, 1}, 'Linear');

Layers = [Layers L];

Controller = FFNN(Layers); % feedforward neural network controller
Plant = DLinearODE(sysd.A, sysd.B, sysd.C, sysd.D, sysd.Ts);
feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial condition of the Plant

% x = [x_lead v_lead a_lead acc2 x_ego v_ego]'

lb = [49; 25; 9; 20];
ub = [51; 25.2; 11; 20.2];

B1 = Box(lb, ub);

init_set = B1.toStar();

%init_set = Polyhedron('lb', lb, 'ub', ub);

% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

lb = [1.4; 30];
ub = [1.4; 31];
B1 = Box(lb, ub);

input_ref = B1.toStar();

%input_ref = Polyhedron('lb', [1.4; 30], 'ub', [1.4; 31]);


%N = 10; % take 8 seconds
%N = 20; % take 16 seconds 
%N = 40; % take 85 seconds, 
N = 40;


n_cores = 4; % number of cores 



% plot output relative distance, ego_car velocity and safe_distance
% distance between two cars
% dis = x_lead - x_ego 

% safe distance between two cars, see here 
% https://www.mathworks.com/help/mpc/examples/design-an-adaptive-cruise-control-system-using-model-predictive-control.html

% dis_safe = D_default + t_gap * v_ego;
% the safety specification is: dis >= 0.9 * dis_safe

% x = [x_lead v_lead a_lead acc2 x_ego v_ego]'
t_gap = 1.4;
D_default = 10;


% safety specification: distance > alp * safe_distance

map = [1 0 -1 0; 0 0 0 1]; % get distance between two cars and velocity of ego car

% plot safe_distance vs. velocity
alp = 1;
map1 = [0 0 0 alp*t_gap; 0 0 0 1]; % safe distance and velocity of ego car
S1 = [];

% figure created before hand as video recording is started, see init file
%figure;
%hold on;
figure(pReach);
hold on;
title('Ego car reduces speed to keep a safe distance with the lead car');
xlabel('Safe distances (red) and actual distances (blue)');
ylabel('Velocity of ego car');
axis equal;
xlim([0 60]);
ylim([-5 25]);

pause(); % wait for user to proceed, so video recording is set up

S = [];

tic(); % start timing

for i=1:N+1
    start_time_i = toc();
    
    if i == 1 % reach(obj, method, init_set, ref_inputSet, n_cores, n_steps)
        [P1, reachTime1] = ncs.reach('approx-star', init_set, input_ref, n_cores, 1);
    else
        [P1, reachTime1] = ncs.reach('approx-star', P1, input_ref, n_cores, 1);
    end
    
    S = [S P1(i).affineMap(map, [])];
    S1 = [S1 P1(i).affineMap(map1, [alp*D_default; 0])];
        
    PLOT.plot_2d_box_noFill(S1, 1, 2, 'red');
    PLOT.plot_2d_box_noFill(S, 1, 2, 'blue');
    end_time_i = toc();
    runtimes(i) = end_time_i - start_time_i;
    ['Time for iteration k = ', num2str(i), ' is: ', num2str(runtimes(i))]
    pause(0.01); % necessary for plot to update
end

['Total runtime for ', num2str(N), ' iterations: ', num2str(toc())]

figure(pRuntime);
hold on;
plot([0:N], runtimes);
xlabel('iteration k');
ylabel('Runtime (s)');
xlim([0 N+1]);
ylim([0 ceil(max(runtimes))*10]);

