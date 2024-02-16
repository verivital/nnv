%% Process some of the results

% for now, we just want to process the multiple attack 3D results
% The results files are 5D 
%   [# images, results, attacktype, # pixels changed, perturbation size]

% res_files = [...
%     "results/verification_multipleAttacks_adrenalmnist3d.mat";
%     "results/verification_multipleAttacks_fracturemnist3d.mat";
%     "results/verification_multipleAttacks_nodulemnist3d.mat";
%     "results/verification_multipleAttacks_organmnist3d.mat";
%     "results/verification_multipleAttacks_synapsemnist3d.mat";
%     "results/verification_multipleAttacks_vesselmnist3d.mat"];

res_files = "results/verification_multipleAttacks_nodulemnist3d.mat";

nFiles = length(res_files);

% Initialize var to summarize results
res_summary = zeros(nFiles, 5);

% Write results to a file
% diary multipleAttacks_3D_results.txt

% Begin processing
for i=1:nFiles

    % Load data
    load(res_files(i));

    % Results title
    disp("*******************************************************");
    disp("================= PROCESSING RESULTS: " + res_files(i) + " ...");

    % Process results for all atacks in file
    process_multiple_attacks(results);

    % Add end of table
    disp("*******************************************************");

end

% Close results file
% diary off;


%% Helper functions

function process_multiple_attacks(results)
    
    % ensure results: 5D array
    n = size(results);
    if length(n) ~= 5
        error("Wrong input.")
    end

    % remember attack combinations
    names = ["dark";"bright"];
    max_pixels = [50;100;200];
    noise_vals = [1;2;3];

    % Begin processing results one attack at a time
    for i = 1:n(3)         % names (dark/bright)
        for j = 1:n(4)     % max_pixels (50/100/200)
            for k = 1:n(5) % noise_val (perturbation size)

                % Print what results we are looking into
                disp("... processing " + names(i) + " attack  with " ...
                    +string(max_pixels(j)) + " pixels perturbed with noise of "+ string(noise_vals(k)));

                % process individual results
                print_results(results(:,:,i,j,k));

                % add empty line for readability
                disp(" ");

            end % k (perturbation size)
        end     % j (# pixels)
    end         % i (attack type)

end

% Print results to screen
function res_summary = print_results(results)
    
    % print results for every attack
    N = size(results,2); 
    disp("----------- ROBUSTNESS RESULTS -------------");
    disp("Verification results of " + string(N) + " images.");
    
    % 1) Robust
    rob = sum(results(1,:) == 1); 
    disp("Number of robust images         =  " + string(rob));
    
    % 0) Not Robust
    not_rob = sum(results(1,:) == 0); 
    disp("Number of not robust images     =  " + string(not_rob));
    
    % 2) Unknown
    unk = sum(results(1,:) == 2);
    disp("Number of unknown images        =  " + string(unk));
    
    % -1) Misclassified
    missclass = sum(results(1,:) == -1); 
    disp("Number of missclassified images =  " + string(missclass));
    
    % Average computation time
    avgTime = sum(results(2,:))/N; 
    disp("Average computation time of        " + string(avgTime));
    
    % Get summary results 
    res_summary = [rob, not_rob, unk, missclass];

end

