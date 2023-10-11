%% Script to verify ACAS Xu properties

% Get path to ACAS Xu data
acas_path = [nnvroot(), filesep, 'data', filesep, 'ACASXu', filesep];

% Iterate through all the networks to verify
networks = dir(acas_path+"onnx/*.onnx");

% property to verify
vnnlibs = dir(acas_path+"vnnlib/*.vnnlib");

% Select network to verify
model_number = 45; % select network 1 out of 45
net_path = [networks(1).folder filesep networks(model_number).name];
net = importNetworkFromONNX(net_path, InputDataFormats='BCSS');

% Select property to verify
prop_name = 'prop_3.vnnlib'; % select property 3 (options 1 to 10)

vnnlib_file = [vnnlibs(1).folder filesep prop_name];

% transform into NNV
net = matlab2nnv(net);

% Define reachability parameters
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

% Verify network
t = tic;
res = net.verify_vnnlib(vnnlib_file, reachOptions);
toc(t);


% Understand model and property
property = load_vnnlib(vnnlib_file); % see help load_vnnlib
% input bounds
input_bounds = [property.lb property.ub]
% property (defined as HalfSpace)
syms x1 x2 x3 x4 x5;
G = property.prop{1}.Hg.G;
g = property.prop{1}.Hg.g;

%  Unsafe if x1 is minimal
property_to_verify = G*[x1;x2;x3;x4;x5] <= g


% How to the output bounds of the network look like?
% Get output set
R = net.reachSet{end}.toStar;
% get ranges of output set
[out_lb, out_ub] = R.getRanges;
% show output bounds
output_bounds = [out_lb, out_ub]

% property is defined as unsafe region
% if property is met, then network is unsafe
% if property is not satisfied, then network is safe

