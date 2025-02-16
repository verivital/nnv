%% Verify all possible 2D classification models for medmnist data

medmnist_path = "data/"; % path to data
models_path = "models/"; % path to trained models

datasets = dir(medmnist_path+"*.mat");

for i=1:length(datasets)

    % get current dataset to verify
    dataset = medmnist_path + datasets(i).name;

    disp("Begin verification of " + datasets(i).name);

    % Load data
    load(dataset);

    % data to verify (test set)
    test_images = permute(test_images, [2 3 4 1]);
    test_labels = test_labels + 1;

    % load network
    load(models_path+"model_"+string(datasets(i).name));
    net = matlab2nnv(net);

    % adversarial attack
    adv_attack = struct;
    adv_attack.Name = "linf";
    adv_attack.epsilon = 1; % {epsilon} color values

    % select images to verify
    % N = 50;
    N = 5;
    inputs = test_images(:,:,:,1:N);
    targets = test_labels(1:N);

    % verify images
    results = verifyDataset(net, inputs, targets, adv_attack);

    % save results
    save("results/verification_"+datasets(i).name, "results", "adv_attack");

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

end