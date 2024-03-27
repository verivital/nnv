%% Verify organ dataset
% VesselMNIST3D, input size: 28x28x28

% Data
dataset = "../../../../../data/medmnist/mat_files/organmnist3d.mat"; % path to data
modelpath = "../../../../../data/medmnist/models/model_organmnist3d.mat";

disp("Begin verification of organ3d");

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
N = 24;
inputs = single(test_images(:,:,:,:,1:N));
targets = single(test_labels(1:N));

% Reachability parameters
reachOptions = struct;
reachOptions.reachMethod = 'relax-star-area';
reachOptions.relaxFactor = 0.95;


%% Attack 1

adv_attack = struct;
adv_attack.Name = "bright";
adv_attack.threshold = 100; % perturb pixels below this value
adv_attack.max_pixels = 50; % Max number of pixels to modify from input image
adv_attack.noise_de = 1; % disturbance (noise) on pixels

results = zeros(N,2); % verification result, time

% verify volumes
parfor i=1:N
    img = squeeze(inputs(:,:,:,:,i));
    target = targets(i);
    results(i,:) = verify_instance_3d(net, img, target, adv_attack, reachOptions);
end

% save results
save("results/verification_organ_"+adv_attack.Name+"_" +adv_attack.noise_de +"_" +adv_attack.max_pixels + ".mat", "results");


%% Attack 2

adv_attack = struct;
adv_attack.Name = "bright";
adv_attack.threshold = 100; % perturb pixels below this value
adv_attack.max_pixels = 100; % Max number of pixels to modify from input image
adv_attack.noise_de = 1; % disturbance (noise) on pixels

results = zeros(N,2); % verification result, time

% verify volumes
parfor i=1:N
    img = squeeze(inputs(:,:,:,:,i));
    target = targets(i);
    results(i,:) = verify_instance_3d(net, img, target, adv_attack, reachOptions);
end

% save results
save("results/verification_organ_"+adv_attack.Name+"_" +adv_attack.noise_de +"_" +adv_attack.max_pixels + ".mat", "results");


%% Attack 3

adv_attack = struct;
adv_attack.Name = "dark";
adv_attack.threshold = 150; % perturb pixels below this value
adv_attack.max_pixels = 50; % Max number of pixels to modify from input image
adv_attack.noise_de = 1; % disturbance (noise) on pixels

results = zeros(N,2); % verification result, time

% verify volumes
parfor i=1:N
    img = squeeze(inputs(:,:,:,:,i));
    target = targets(i);
    results(i,:) = verify_instance_3d(net, img, target, adv_attack, reachOptions);
end

% save results
save("results/verification_organ_"+adv_attack.Name+"_" +adv_attack.noise_de +"_" +adv_attack.max_pixels + ".mat", "results");


%% Attack 4

adv_attack = struct;
adv_attack.Name = "dark";
adv_attack.threshold = 150; % perturb pixels below this value
adv_attack.max_pixels = 100; % Max number of pixels to modify from input image
adv_attack.noise_de = 1; % disturbance (noise) on pixels

results = zeros(N,2); % verification result, time

% verify volumes
parfor i=1:N
    img = squeeze(inputs(:,:,:,:,i));
    target = targets(i);
    results(i,:) = verify_instance_3d(net, img, target, adv_attack, reachOptions);
end

% save results
save("results/verification_organ_"+adv_attack.Name+"_" +adv_attack.noise_de +"_" +adv_attack.max_pixels + ".mat", "results");


%% Attack 5

adv_attack = struct;
adv_attack.Name = "bright";
adv_attack.threshold = 100; % perturb pixels below this value
adv_attack.max_pixels = 50; % Max number of pixels to modify from input image
adv_attack.noise_de = 2; % disturbance (noise) on pixels

results = zeros(N,2); % verification result, time

% verify volumes
parfor i=1:N
    img = squeeze(inputs(:,:,:,:,i));
    target = targets(i);
    results(i,:) = verify_instance_3d(net, img, target, adv_attack, reachOptions);
end

% save results
save("results/verification_organ_"+adv_attack.Name+"_" +adv_attack.noise_de +"_" +adv_attack.max_pixels + ".mat", "results");


%% Attack 6

adv_attack = struct;
adv_attack.Name = "bright";
adv_attack.threshold = 100; % perturb pixels below this value
adv_attack.max_pixels = 100; % Max number of pixels to modify from input image
adv_attack.noise_de = 2; % disturbance (noise) on pixels

results = zeros(N,2); % verification result, time

% verify volumes
parfor i=1:N
    img = squeeze(inputs(:,:,:,:,i));
    target = targets(i);
    results(i,:) = verify_instance_3d(net, img, target, adv_attack, reachOptions);
end

% save results
save("results/verification_organ_"+adv_attack.Name+"_" +adv_attack.noise_de +"_" +adv_attack.max_pixels + ".mat", "results");


%% Attack 7

adv_attack = struct;
adv_attack.Name = "dark";
adv_attack.threshold = 150; % perturb pixels below this value
adv_attack.max_pixels = 50; % Max number of pixels to modify from input image
adv_attack.noise_de = 2; % disturbance (noise) on pixels

results = zeros(N,2); % verification result, time

% verify volumes
parfor i=1:N
    img = squeeze(inputs(:,:,:,:,i));
    target = targets(i);
    results(i,:) = verify_instance_3d(net, img, target, adv_attack, reachOptions);
end

% save results
save("results/verification_organ_"+adv_attack.Name+"_" +adv_attack.noise_de +"_" +adv_attack.max_pixels + ".mat", "results");


%% Attack 8

adv_attack = struct;
adv_attack.Name = "dark";
adv_attack.threshold = 150; % perturb pixels below this value
adv_attack.max_pixels = 100; % Max number of pixels to modify from input image
adv_attack.noise_de = 2; % disturbance (noise) on pixels

results = zeros(N,2); % verification result, time

% verify volumes
parfor i=1:N
    img = squeeze(inputs(:,:,:,:,i));
    target = targets(i);
    results(i,:) = verify_instance_3d(net, img, target, adv_attack, reachOptions);
end

% save results
save("results/verification_organ_"+adv_attack.Name+"_" +adv_attack.noise_de +"_" +adv_attack.max_pixels + ".mat", "results");

