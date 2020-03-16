% /* An example of verifying a continuous linear NNCS */
% /* Controller
load controller_5_20.mat; weights = network.weights;
bias = network.bias; n = length(weights); Layers = [];
for i=1:n - 1
    L = LayerS(weights{1, i}, bias{i, 1}, 'poslin');
    Layers = [Layers L];
end
L = LayerS(weights{1, n}, bias{n, 1}, 'purelin');
Layers = [Layers L];
Controller = FFNNS(Layers); 
% /* plant model
A = [0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 0 0 0 1; ...
    0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 -2 0; ...
    0 0 -2 0 0 0 0];
B = [0; 0; 0; 0; 0; 2; 0];
C = [1 0 0 -1 0 0 0; 0 1 0 0 -1 0 0; 0 0 0 0 1 0 0];  
D = [0; 0; 0];
Tc = 0.1; % control period
Nr = 20; % number of reachability steps in 1 control period
plant = LinearODE(A, B, C, D, Tc, Nr); % continuous plant model
% /* continuous linear NNCS 
ncs = LinearNNCS(Controller, plant); % a continuous linear NNCS
% /* ranges of initial set of states of the plant
lb = [90; 29; 0; 30; 30; 0; -10];
ub = [92; 30; 0; 31; 30.5; 0; -10];
% /* reachability parameters
reachPRM.init_set = Star(lb, ub);
reachPRM.ref_input = [30; 1.4];
reachPRM.numSteps = 10;
reachPRM.reachMethod = 'approx-star';
reachPRM.numCores = 4;
% /* usafe region: x1 - x4 <= 1.4 * v_ego + 10
unsafe_mat = [1 0 0 -1 -1.4 0 0];
unsafe_vec = 10;
U = HalfSpace(unsafe_mat, unsafe_vec);
%U = HalfSpace([-1 0 0 0 0 0 0], -20);

% /* verify the system
[safe, counterExamples, verifyTime] = ncs.verify(reachPRM, U);

