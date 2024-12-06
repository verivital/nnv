% Run smoke test.

datasets = ["zoom_in", "zoom_out", "gtsrb", "stminst"];
veralgs = ["relax", "approx"]
timeout = 1800
epsilon_index = 1;

% Compute number of samples to verify in smoke test
numSamples = length(datasets) * length(veralgs) * 3;

% Make a matrix to save results
results = zeros([numSamples+1 3]);

% Add headers to the results file
results(1, 1) = "Res";
results(1, 2) = "Time";
results(1, 3) = "Met";

num_samples_verified = 2;

for dataset_index=1:length(datasets)
    dataset = datasets(dataset_index);

    if dataset == "stmnist"
        frame_nums = [16, 32, 64];
    else
        frame_nums = [4, 8, 16];
    end

    for frame_num_index=1:length(frame_nums)
        frame_num = frame_nums(frame_num_index);

        for veralg_index=1:length(veralgs)
            veralg = veralgs(veralg_index);

            if dataset == "zoom_in" || dataset == "zoom_out"
                if frame_num == 16
                    continue
                end
                [res, t, met] = verifyvideo(dataset, frame_num, veralg, 1, 751)

            elseif dataset == "gtsrb"
                [res, t, met] = verifygtsrb(dataset, frame_num, veralg, 1, 1)

            else
                [res, t, met] = verifystmnist(dataset, frame_num, veralg, 1, 1)
            end
            
            % Add the results to the table
            results(num_samples_verified, 1) = res;
            results(num_samples_verified, 2) = t;
            results(num_samples_verified, 3) = met;

            % Increment to next row in results table
            num_samples_verified = num_samples_verified + 1;
        end
    end
end

% Save the results
filename = "smoketest_outputs.txt";
writematrix(results, filename);
