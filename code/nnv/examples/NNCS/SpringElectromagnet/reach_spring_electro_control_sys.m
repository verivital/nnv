load controller.mat;
load sysd.mat;

L1 = LayerS(W{1,1}, b{1,1}, 'poslin');
L2 = LayerS(W{1,2}, b{1,2}, 'purelin');

NN_Controller = FFNNS([L1 L2]); % feedforward neural network controller
Plant = DLinearODE(sysd.A, sysd.B, sysd.C, sysd.D, 0.1);
feedbackMap = [0;1]; % feedback map, y[k] and y[k-1]

ncs = NNCS(NN_Controller, Plant, feedbackMap); % the neural network control system

init_set = Star([0.8;-1],[1;-0.8]);
input_set = Star(0.99,1);
% input_set = 1;
% % initial set of state of the plant 
% init_set = Polyhedron('lb', [0.8; -1], 'ub', [1; -0.8]);
% 
% % control input set: 0.99 m <= r <= 1 m
% ref_inputSet = Polyhedron('lb', 0.99, 'ub', 1);
% 
N = 5; % number of step
n_cores = 1; % number of cores 
reachPRM.numCores = n_cores;
reachPRM.init_set = init_set;
reachPRM.numSteps = N;
reachPRM.reachMethod = 'approx-star';
reachPRM.ref_input = input_set;
reachPRM2 = reachPRM;
reachPRM2.reachMethod = 'exact-star';
[P1, reachTime1] = ncs.reach(reachPRM);
[P2, reachTime2] = ncs.reach(reachPRM2);

% 
% [P1, reachTime1] = ncs.reach('approx-polytope', init_set, ref_inputSet, n_cores, N);
% [P2, reachTime2] = ncs.reach('exact-polytope', init_set, ref_inputSet, n_cores, N);
% 
% figure;
% P1.plot;
% figure;
% P2.plot;