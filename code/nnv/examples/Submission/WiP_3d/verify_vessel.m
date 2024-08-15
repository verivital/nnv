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

% Study variables
advType = ["add", "remove"];
maxpixels = [50, 100, 500, 1000]; %out of 28x28x28 pixels
epsilon = 1; % ep / 255


%% Verification analysis
for a=advType
    for mp=maxpixels
        for ep=epsilon
            
            % 1) Initialize results var
            results = zeros(N,2);
            
            % 2) Create adversarial attack
            adv_attack = struct;
            adv_attack.Name = a; % add or remove
            adv_attack.max_pixels = mp; % Max number of pixels to modify from input image
            adv_attack.noise_de = ep/255; % disturbance (noise) on pixels
            
            % 3) Begin verification analysis
            for i=1:N
                img = squeeze(inputs(:,:,:,:,i));
                target = targets(i);
                results(i,:) = verify_instance_shape(net, img, target, adv_attack, reachOptions);
            end
            
            % 4) % save results
            save("results/verification_vessel_"+ a +"_" + ep +"_" + mp + ".mat", "results");

        end
    end
end

