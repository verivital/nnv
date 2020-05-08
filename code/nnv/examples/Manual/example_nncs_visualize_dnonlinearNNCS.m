% /* An example of constructing a discrete nonlinear NNCS */
%/* FFNN controller
load MountainCar_ReluController.mat;
W = nnetwork.W; % weight matrices
b = nnetwork.b; % bias vectors
n = length(W);
Layers = [];
for i=1:n - 1
    L = LayerS(W{1, i}, b{1, i}, 'poslin');
    Layers = [Layers L];
end
L = LayerS(W{1, n}, b{1, n}, 'purelin');
Layers = [Layers L];
controller = FFNNS(Layers);

%/* MountainCar
Ts = 0.5; % sampling time
C = [1 0; 0 1]; % output matrix
Car = DNonLinearODE(2, 1, @discrete_car_dynamics, Ts, C); 

ncs = DNonlinearNNCS(controller, Car); % system

lb = [-0.41; 0];
ub = [0.4; 0];
reachPRM.init_set = Star(lb, ub);
reachPRM.ref_input = [];
reachPRM.numSteps = 10;
reachPRM.reachMethod = 'approx-star';
reachPRM.numCores = 4;

% unsafe region
U = HalfSpace([-1 0], 0); % x1 > 0

[safe, counterExamples, verifyTime] = ncs.verify(reachPRM, U);

%% plot 2d output sets
figure; 
map_mat = [1 0; 0 1];
map_vec = [];
ncs.plotOutputReachSets('blue', map_mat, map_vec);
title('Position versus speed');

%% plot 1d output sets
figure; 
map_mat = [1 0];
map_vec = [];
ncs.plotOutputReachSets('blue', map_mat, map_vec);
title('Position versus time');
