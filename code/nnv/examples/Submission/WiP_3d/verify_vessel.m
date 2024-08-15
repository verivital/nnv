%% Verify vessel dataset
% vesselMNIST3D, input size: 28x28x28

% Data
dataset = "../../../../../data/medmnist/mat_files/vesselmnist3d.mat"; % path to data
modelpath = "../../../../../data/medmnist/models/model_vesselmnist3d.mat";

disp("Begin verification of vessel3d");

% Load data
load(dataset);

% data to verify (test set)
test_images = permute(test_images, [2 3 4 5 1]);
test_labels = test_labels + 1;

% load network
load(modelpath);
matlabNet = net;
net = matlab2nnv(net);

% select volumes to verify
N = 24; % even for numCores
inputs = single(test_images(:,:,:,:,1:N));
targets = single(test_labels(1:N));

% Reachability parameters
reachOptions = struct;
reachOptions.reachMethod = 'relax-star-area';
reachOptions.relaxFactor = 0.95;
reachOptions.lp_solver = "gurobi";

%% Attack 1

% adv_attack = struct;
%???????

% results = zeros(N,2); % verification result, time

