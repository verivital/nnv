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
N = 40; % take 85 seconds, 


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

map = [1 0 -1 0; 0 0 0 1]; % get distance between two cars and velocity of ego car

S = [];
for i=1:N+1
    S = [S P1(i).affineMap(map, [])];
end


% plot safe_distance vs. velocity
alp = 1;
map1 = [0 0 0 alp*t_gap; 0 0 0 1]; % safe distance and velocity of ego car
S1 = [];
for i=1:N+1
    S1 = [S1 P1(i).affineMap(map1, [alp*D_default; 0])];
end

figure;
PLOT.plot_2d_box_noFill(S, 1, 2, 'blue');
title('Ego car reduces speed to keep a safe distance with the lead car');
xlabel('Safe distances (red) and actual distances (blue)');
ylabel('Velocity of ego car');

hold on;
PLOT.plot_2d_box_noFill(S1, 1, 2, 'red');




