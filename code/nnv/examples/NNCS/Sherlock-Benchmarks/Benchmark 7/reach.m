% Notes
% This examples errors out (CORA error in initReach)
%
%
% Index in position 2 exceeds array bounds. Index must not exceed 2.
% 
% Error in zonotope/split>splitOneDim (line 106)
% c1=c-G(:,dim)/2;
% 
% Error in zonotope/split (line 49)
%         Zsplit{N}=splitOneDim(Z,N);
% 
% Error in select (line 34)
% Rtmp = split(Rinit.set);
% 
% Error in contDynamics/linReach (line 180)
%     dimForSplit = select(obj,options,Rstart);
% 
% Error in nonlinearSys/initReach (line 52)
%         [Rtemp_ti,Rtemp_tp,dimForSplit,opts] = linReach(obj,options,Rinit{i});
%
% Error in nonlinearSys/initReach (line 80)
%             Rres = initReach(obj,Rsplit,options);
% 
% Error in contDynamics/reach (line 120)
%     [Rnext, options] = initReach(obj,options.R0,options);
% 
% Error in NonLinearODE/reach_zono (line 240)
%            R = reach(sys, obj.params, obj.options); % CORA reach method using zonotope and conservative linearization
% 
% Error in NonLinearODE/stepReachStar (line 327)
%                 [R, ~] = obj.reach_zono(I, U, obj.options.timeStep, obj.params.tFinal);
% 
% Error in NNCS/reach (line 169)
%                  R = obj.plant.stepReachStar(fb_I{1}(length(fb_I{1})), U1);
% 
% Error in reach (line 50)
% [P, reachTime] = ncs.reach(reachPRM);

%% NNCS

% Controller
load controller.mat;

weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = {};
for i=1:n-1
    Layers{i} = LayerS(weights{1, i}, bias{1, i}, 'poslin');
end
Layers{n} = LayerS(weights{1, n}, bias{1, n}, 'purelin');
Controller = NN(Layers); % feedforward neural network controller
Controller.InputSize = 3;
Controller.OutputSize = 1;

%Plant
reachStep = 0.02;
controlPeriod = 0.5;
output_mat = eye(3); % feedback 
Plant = NonLinearODE(3, 1, @dynamics, reachStep, controlPeriod, output_mat);

% NNCS
feedbackMap = [0]; % feedback map, y[k] 
ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system


%% Analysis
% N = 50; % number of control steps   
N = 3;

% lb = [0.35; 0.45; 0.25];
lb = [0.445; 0.50; 0.30];
ub = [0.45; 0.55; 0.35];
init_set = Star(lb, ub);
input_ref = [];

% Evaliation
[simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, []);

% Reachability parameters
reachPRM.numSteps = N;
reachPRM.numCores = 1;
reachPRM.reachMethod = 'approx-star';
reachPRM.ref_input = input_ref;
reachPRM.init_set = init_set;

% Reachability computation
[P, reachTime] = ncs.reach(reachPRM);


%% Visualization
fig = figure;
Star.plotBoxes_2D_noFill(P, 1, 2, 'blue');


%% Notes