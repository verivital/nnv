%% Debug inconsistencies for nodulemnist verification results

% get current dataset to verify
% dataset = medmnist_path + datasets(i).name;
dataset = "data/mat_files/nodulemnist3d.mat";


disp("Begin verification of nodulemnist3d");

% Load data
load(dataset);

% data to verify (test set)
test_images = permute(test_images, [2 3 4 5 1]);
test_labels = test_labels + 1;

% load network
load("models/model_nodulemnist3d.mat");
matlabNet = net;
net = matlab2nnv(net);

% adversarial attacks
names = ["dark"];
max_pixels = [50;100;200];
noise_vals = [1];

% select volumes to verify
N = 200;
inputs = test_images(:,:,:,:,1:N);
targets = test_labels(1:N);

% Initialize results 
results = zeros(2, N, length(names), length(max_pixels), length(noise_vals));
outputSets = cell(3,1);

% verify volumes with all attack combos
for a=1:length(names)
    for b=1:length(max_pixels)
        for c=1:length(noise_vals)
            % create attack from variables
            adv_attack = struct;
            adv_attack.Name = names(a);
            adv_attack.max_pixels = max_pixels(b);
            adv_attack.noise_de = noise_vals(c);
            adv_attack.threshold = 150;

            % Compute verification
            [outputSets{b}, results(:,:,a,b,c)] = verify_medmnist3d_extraInfo(net, inputs, targets, adv_attack);
        end
    end
end

% save results
save("results/verification_multipleAttacks_nodulemnist3d_debug.mat", "results", "outputSets");


