%% Verify all possible 3D classification models for medmnist data


medmnist_path = "data/mat_files/"; % path to data

datasets = dir(medmnist_path+"*.mat");

for i=1:length(datasets)

    if endsWith(datasets(i).name, "3d.mat")

        % get current dataset to verify
        dataset = medmnist_path + datasets(i).name;

        try % verification

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

            % adversarial attack
            attack_option = 'bright'; % default (darkening attack)
            switch attack_option
                case 'linf'
                    adv_attack = struct;
                    adv_attack.Name = "linf";
                    adv_attack.epsilon = 1; % {epsilon} color values
                    % adv_attack.max_pixels = 784; % Max number of pixels to modify from input image
                    adv_attack.max_pixels = 1; % Max number of pixels to modify from input image
                case 'dark'
                    adv_attack = struct;
                    adv_attack.Name = "dark";
                    adv_attack.threshold = 150; % perturb pixels over this value
                    adv_attack.max_pixels = 50; % Max number of pixels to modify from input image
                    adv_attack.noise_de = 1; % disturbance (noise) on pixels
                case 'bright'
                    adv_attack = struct;
                    adv_attack.Name = "bright";
                    adv_attack.threshold = 100; % perturb pixels below this value
                    adv_attack.max_pixels = 50; % Max number of pixels to modify from input image
                    adv_attack.noise_de = 1; % disturbance (noise) on pixels
                otherwise
                    error("Wrong attack");
            end

            % select volumes to verify
            N = 200;
            inputs = test_images(:,:,:,:,1:N);
            targets = test_labels(1:N);

            % verify volumes
            results = verify_medmnist3d(net, inputs, targets, adv_attack);

            % save results
            save("results/verification_"+adv_attack.Name+"_"+datasets(i).name, "results", "adv_attack");

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

        catch ME

            warning("Failed!!")
            warning(ME.message);
            disp(medmnist_path+datasets(i).name)
        end

    end

end