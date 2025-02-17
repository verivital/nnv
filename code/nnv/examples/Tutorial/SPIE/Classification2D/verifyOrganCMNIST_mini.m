%% Verify all possible 2D classification models for medmnist data

medmnist_path = "data/organcmnist.mat"; % path to data
model_path = "models/model_organcmnist.mat"; % path to trained models

disp("Begin verification of organCMNIST");

% Load data
load(medmnist_path);

% data to verify (test set)
test_images = permute(test_images, [2 3 4 1]);
test_labels = test_labels + 1;

% load network
load(model_path);
net = matlab2nnv(net);

% select images to verify
% N = 50;
N = 10;
inputs = test_images(:,:,:,1:N);
targets = test_labels(1:N);

% adversarial attack
adv_attack = struct;
adv_attack.Name = "linf";
adv_attack.epsilon = 1; % {epsilon} color values

% verify images
results = verifyDataset(net, inputs, targets, adv_attack);

% save results
save("results/verification_organcmnist_mini.mat", "results", "adv_attack");

% print results to screen
disp("======= ROBUSTNESS RESULTS ==========")
disp(" ");
disp("Verification results of " + string(N) + " images.")
disp("Number of robust images          =  " + string(sum(results(1,:) == 1)));
disp("Number of not robust images      =  " + string(sum(results(1,:) == 0)));
disp("Number of unknown images         =  " + string(sum(results(1,:) == 2)));
disp("Number of missclassified images  =  " + string(sum(results(1,:) == -1)))
disp(" ");
disp("Total computation time of " + string(sum(results(2,:))));

disp("|========================================================================|")
disp(' ');

