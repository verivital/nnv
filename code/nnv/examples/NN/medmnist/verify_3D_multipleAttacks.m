%% Verify all possible 3D classification models for medmnist data

% Path to datasets
medmnist_path = "data/mat_files/"; % path to data
datasets = dir(medmnist_path+"*3d.mat");

% Define reachability parameters
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
% reachOptions.reachMethod = 'relax-star-area';
% reachOptions.relaxFactor = 0.5; % solve only 1-relaxFactor of LPs for the 2 relu layers in the network

for i=1:length(datasets)

    % get current dataset to verify
    dataset = medmnist_path + datasets(i).name;

    disp("Begin verification of " + datasets(i).name);

    % Load data
    load(dataset);

    % data to verify (test set)
    test_images = permute(test_images, [2 3 4 5 1]);
    test_labels = test_labels + 1;

    % load network
    load("models/model_"+string(datasets(i).name));
    matlabNet = net;
    net = matlab2nnv(net);

    % adversarial attacks
    names = ["dark", "bright","linf"];
    max_pixels = [50; 100; 200];
    noise_vals = [1/255; 2/255; 3/255];

    % select volumes to verify
    N = 200;
    inputs = test_images(:,:,:,:,1:N);
    targets = test_labels(1:N);

    % Initialize results 
    results = zeros(2, N, length(names), length(max_pixels), length(noise_vals));
    errors = cell(length(names), length(max_pixels), length(noise_vals));

    % verify volumes with all attack combos
    for a = 1:length(names)
        for b = 1:length(max_pixels)
            for c = 1:length(noise_vals)
                % create attack from variables
                adv_attack = struct;
                adv_attack.Name = names(a);
                adv_attack.max_pixels = max_pixels(b);
                adv_attack.noise_de = noise_vals(c);
                if strcmp(names(a),'dark')
                    adv_attack.threshold = 150;
                elseif strcmp(names(a), 'bright')
                    adv_attack.threshold = 100;
                end
                fprintf("Reach set of %s with %d pixels, threshold %d , and noise %f \n",...
                    adv_attack.Name, adv_attack.max_pixels, adv_attack.threshold, adv_attack.noise_de);
                % Compute verification
                [results(:,:,a,b,c), errors{a,b,c}] = verify_medmnist3d(net, inputs, targets, adv_attack, reachOptions);
            end
        end
    end

    % save results
    save("results/verification_darkening_"+datasets(i).name, "results", "reachOptions");

end

