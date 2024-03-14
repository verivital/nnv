%% Fairness Verification of Adult Classification Model (NN)

% Load network 
ac_model = load('adult_model_fc.mat');

% Create NNV model 
net = matlab2nnv(ac_model.net);

%% Verification
% Utilize the verify_robustness function?
% Need to initalize IS (input set), reachOptions, and target

% First, we define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approxiate method

% Second, we define the input set


%Third, we define the target



