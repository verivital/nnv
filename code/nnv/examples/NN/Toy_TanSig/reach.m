% Load network
load('Engine_Toy_Tansig_net.mat');

% Get parameters from network
% bias
b = net.b;
nL = length(b); % number of layers
% weights
W = cell(nL,1);
W{1} = net.IW{1};
W{2} = net.LW{2,1};
W{3} = net.LW{3,2};
% Create NN model for NNV
Layers = cell(nL,1);
for i=1:nL
    Layers{i} = LayerS(W{i}, b{i}, net.Layers{i}.transferFcn);
end
% NN object for the reachability analysis
nn = NN(Layers); 

% Define input set
lb = [-1; -1];
ub = [1; 1];
I = Star(lb, ub); % create input set as a Star set

% Compute reachability
fprintf('\nComputing the reachable set of the toy network:\n');
reachOptions = struct; % define reachability options
reachOptions.reachMethod = 'approx-star'; % define reachability method
reachSet = nn.reach(I, reachOptions); % perform reachability analyis using approx-star
[y,Y] = reachSet.getRanges; % get the output ranges (overapproximation)

disp('Example completed');
