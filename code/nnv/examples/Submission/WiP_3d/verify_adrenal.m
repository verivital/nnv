%% Verify adrenal dataset
% adrenalMNIST3D, input size: 28x28x28

% Data
dataset = "../../../../../data/medmnist/mat_files/adrenalmnist3d.mat"; % path to data
modelpath = "../../../../../data/medmnist/models/model_adrenalmnist3d.mat";

disp("Begin verification of adrenal3d");

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
% N = 24; % even for numCores
N = 50;
inputs = single(test_images(:,:,:,:,1:N));
targets = single(test_labels(1:N));

% Reachability parameters
reachOptions = struct;
reachOptions.reachMethod = 'relax-star-area';
reachOptions.relaxFactor = 0.95;
reachOptions.lp_solver = "gurobi";


%% Verification analysis
