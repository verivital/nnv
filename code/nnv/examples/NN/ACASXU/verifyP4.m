%% Script to verify ACAS Xu property 4 (all 45 networks)

% Get path to ACAS Xu data
acas_path = [nnvroot(), filesep, 'data', filesep, 'ACASXu', filesep];

% Iterate through all the networks to verify
networks = dir(acas_path+"onnx/*.onnx");

% property to verify
vnnlib_file = acas_path+"vnnlib/prop_4.vnnlib";

% Define reachability parameters
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

% Preallocate memory for results
results = zeros(45,1) - 1;
reachTime = zeros(45,1);

% Begin reachability
for i = 1:length(networks)

    % Load Network
    file = [networks(i).folder, filesep,  networks(i).name];
    net = importONNXNetwork(file, InputDataFormats='BCSS');

    % transform into NNV
    net = matlab2nnv(net);
    
    % Verify network
    Rstar = []; % to transform the ImageStar into Starset (vnnlib properties get automatically defined as ImageStars)
    t = tic;
    results(i) = net.verify_vnnlib(vnnlib_file, reachOptions);
    reachTime(i,1) = toc(t);
    
end
