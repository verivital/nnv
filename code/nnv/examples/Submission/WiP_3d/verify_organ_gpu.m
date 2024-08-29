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
N = 50;
idxs = zeros(N,1);
count = 1;

for i = 1:length(test_labels)
    y = classify(matlabNet,test_images(:,:,:,i));
    if single(y) == test_labels(i)
        idxs(count) = i;
        count = count + 1;
    end
    if count > N
        break
    end
end

inputs = single(test_images(:,:,:,idxs));
targets = single(test_labels(idxs));

% Reachability parameters
reachOptions = struct;
reachOptions.reachMethod = 'relax-star-area';
reachOptions.relaxFactor = 0.95;
reachOptions.lp_solver = "gurobi";
reachOptions.device = 'gpu';

% Study variables
advType = ["bright", "dark"];
maxpixels = [50, 100, 500, 1000]; %out of 28x28x28 pixels
epsilon = [2, 5, 10, 25]; % ep / 255
threshold = [100; 150]; % bright ; dark

%% Verification analysis
for a=advType
    for mp=maxpixels
        for ep=epsilon
            
            % 1) Initialize results var
            results = zeros(N,2);
            
            % 2) Create adversarial attack
            adv_attack = struct;
            adv_attack.Name = a; % bright or dark
            if strcmp(a, "bright") 
                adv_attack.threshold = threshold(1); % perturb pixels below this value
            else 
                adv_attack.threshold = threshold(2); % perturb pixels below this value
            end 
            adv_attack.max_pixels = mp; % Max number of pixels to modify from input image
            adv_attack.noise_de = ep/255; % disturbance (noise) on pixels
            
            % 3) Begin verification analysis
            for i=1:N
                img = inputs(:,:,:,i);
                target = targets(i);
                results(i,:) = verify_instance_3d(net, img, target, adv_attack, reachOptions);
            end
            
            % 4) % save results
            save("resultsGPU/verification_organ_"+ a +"_" + ep +"_" + mp + ".mat", "results");

        end
    end
end