load controller_10_20.mat;
weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = [];
for i=1:n - 1
    L = Layer(weights{1, i}, bias{i, 1}, 'ReLU');
    Layers = [Layers L];
end

L = Layer(weights{1, n}, bias{n, 1}, 'Linear');

Layers = [Layers L];

Controller = FFNN(Layers); % feedforward neural network controller
Plant = NonLinearODE(6, 1, @car_dynamics);
Plant.set_timeStep(0.01); % time step for reachability analysis of the plant
Plant.set_tFinal(0.1); % Ts = 0.1, sampling time for control signal from neural network controller
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial condition of the Plant

% x = [x_lead v_lead x_internal_lead x_ego v_ego x_internal_ego]'

lb = [98; 32; 0; 10; 30; 0];
ub = [100; 32.2; 0; 10; 30.2; 0];

init_set = Star(lb, ub);

% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

lb = [30; 1.4];
ub = [30; 1.4];
B1 = Box(lb, ub);

input_ref = B1.toStar();

N = 50;  

n_cores = 4; % number of cores 

[P1, reachTime1] = ncs.reach('approx-star', init_set, input_ref, n_cores, N);

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

map = [1 0 0 -1 0 0; 0 0 0 0 1 0]; % get distance between two cars and velocity of ego car

S = [];
for i=1:N+1
    S = [S P1(i).affineMap(map, [])];
end


% plot safe_distance vs. velocity
alp = 1;
map1 = [0 0 0 0 alp*t_gap 0; 0 0 0 0 1 0]; % safe distance and velocity of ego car
S1 = [];
for i=1:N+1
    S1 = [S1 P1(i).affineMap(map1, [alp*D_default; 0])];
end

%figure;
%PLOT.plot_2d_box_noFill(S, 1, 2, 'blue'); % plot distance vs velocity
%hold on;
%PLOT.plot_2d_box_noFill(S1, 1, 2, 'red'); % plot safe distance vs. velocity
%title('Ego car reduces speed to keep a safe distance with the lead car');
%xlabel('Safe distances (red) and actual distances (blue)');
%ylabel('Velocity of ego car');


% plot intermediate reach set (all reachable set of the plant) vs. time
reachSet = ncs.plant.intermediate_reachSet;
reachSet = [init_set reachSet]; % add init_set into the reachable set
times = 0:0.01:N*0.1;

dis = [];
safe_dis = [];
ego_vel = [];
lead_vel = [];

for i=1:length(reachSet)
    dis = [dis reachSet(i).affineMap([1 0 0 -1 0 0], [])];
    safe_dis = [safe_dis reachSet(i).affineMap([0 0 0 0 alp*t_gap 0], alp*D_default)];
    ego_vel = [ego_vel reachSet(i).affineMap([0 0 0 0 1 0], [])];
    lead_vel = [lead_vel reachSet(i).affineMap([0 1 0 0 0 0], [])];
end


% get controller output set and plot it
U = ncs.controlSet; 


% plot velocity, distance, safe_distance versus time
figure;
subplot(3, 1, 1);
Star.plotRanges_2D(ego_vel, 1, times, 'blue'); % plot ego car's velocity versus time
hold on;
Star.plotRanges_2D(lead_vel, 1, times, 'green'); % plot ego car's velocity versus time
xlabel('time (seconds)');
ylabel('Velocity');
title('Ego car velocity (blue) vs. lead car velocity (green)');
subplot(3, 1, 2);
Star.plotRanges_2D(dis, 1, times, 'blue'); % plot distance between two cars versus time 
hold on;
Star.plotRanges_2D(safe_dis, 1, times, 'red'); % plot safe distance versus time
xlabel('time(seconds)');
ylabel('Distance');
title('Actual distance (blue) vs. safe distance (red)');
subplot(3, 1, 3);
times1 = 0:0.1:N*0.1;
Star.plotRanges_2D(U, 1, times1(1:N), 'cyan'); % plot control set versus time 
xlabel('time(seconds)');
ylabel('Control');
title('Control vs. time');
