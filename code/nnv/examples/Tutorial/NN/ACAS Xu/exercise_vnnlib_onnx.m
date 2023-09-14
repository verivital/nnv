%% Script to verify ACAS Xu properties

% Get path to ACAS Xu data
acas_path = [nnvroot(), filesep, 'data', filesep, 'ACASXu', filesep];

% Iterate through all the networks to verify
networks = dir(acas_path+"onnx/*.onnx");

% property to verify
vnnlibs = dir(acas_path+"vnnlib/*.vnnlib");

% Select network to verify
model_number = ; % select network (options 1 to 45)
net_path = [networks(model_number).folder filesep networks(model_number).name];
net = importONNXNetwork(); % can use help / doc importONNXNetwork for more information 
% transform into NNV
net = ; % transform MATLAB's network to NNV format

% Select property to verify
prop_name = 'prop_(enter number here).vnnlib'; % select property (options 1, 2, 3, 4 or 5)
vnnlib_file = [vnnlibs(1).folder filesep prop_name];

% Define reachability parameters
reachOptions = struct;
reachOptions.reachMethod = ; % sound and complete (exact) or sound and incomplete (approximate)

% Verify network
t = tic;
res = net.verify_vnnlib(vnnlib_file, reachOptions);
toc(t);
