%% Process some of the results

% for now, we just want to process the multiple attack 3D results
% The results files are 5D ([#images, results, 

res_files = [...
    "results/verification_multipleAttacks_adrenalmnist3d.mat";
    "results/verification_multipleAttacks_fracturemnist3d.mat";
    "results/verification_multipleAttacks_nodulemnist3d.mat";
    "results/verification_multipleAttacks_organmnist3d.mat";
    "results/verification_multipleAttacks_synapsemnist3d.mat";
    "results/verification_multipleAttacks_vesselmnist3d.mat"];

nFiles = length(res_files);
res_summary = zeros(nFiles, 5);
% Write results to a file
diary multipleAttacks_3D_results.txt
for i=1:nFiles
    load(res_files(i));
    disp("*******************************************************");
    disp("================= PROCESSING RESULTS: " + res_files(i) + " ...");
    process_multiple_attacks(results);
    disp("*******************************************************");
end
diary off;


%% Helper functions

function process_multiple_attacks(results)
    % results: 5D array
    n = size(results);
    if length(n) ~= 5
        error("Wrong input.")
    end
    % attack combinations
    names = ["dark";"bright"];
    max_pixels = [50;100;200];
    noise_vals = [1;2;3];
    for i = 1:n(3)
        for j = 1:n(4)
            for k = n(5)
                % Print what results we are looking into
                disp("... processing " + names(i) + " attack  with " ...
                    +string(max_pixels(j)) + " pixels perturbed with noise of "+ string(noise_vals(k)));
                % process individual results
                print_results(results(:,:,i,j,k));
                disp(" ");
            end
        end
    end
end

% Print results to screen
function res_summary = print_results(results)
    % print results for every attack
    N = size(results,2); 
    disp("----------- ROBUSTNESS RESULTS -------------");
    disp("Verification results of " + string(N) + " images.");
    rob = sum(results(1,:) == 1); % robust
    disp("Number of robust images          =  " + string(rob));
    not_rob = sum(results(1,:) == 0); % not robust
    disp("Number of not robust images      =  " + string(not_rob));
    unk = sum(results(1,:) == 2); % unknown
    disp("Number of unknown images         =  " + string(unk));
    missclass = sum(results(1,:) == -1); % misclassified
    disp("Number of missclassified images  =  " + string(missclass));
    avgTime = sum(results(2,:))/N; % average computation time
    disp("Average computation time of "         + string(avgTime));
    % Get summary results 
    res_summary = [rob, not_rob, unk, missclass];
end

