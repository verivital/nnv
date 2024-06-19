%% Get all results, one dataset at a time

results = dir('results');

% Variable in the study
datasets = ["adrenal"; "fracture"; "nodule"; "organ"; "synapse"; "vessel"];
attackType = ["bright"; "dark"];
maxpixels = [50;100;500;1000];
epsilon = [2/255; 4/255; 10/255];

%% Visualize verification results/trends

% What do we want to show?

% 1) Trends
%  - Time: How does the time increases as the input set gets larger?
%     - Epsilon increasing
%     - Number of pixels increasing
%  - results: How robust are the models as we increase the size of the attack? more unknowns? more sats?
%     - Epsilon increasing
%     - Number of pixels increasing
% 2) Certified robust accuracy
%  - Is it much worse than the accuracy of the model? 
%  - Do we need more samples?
 


%% Helper functions

function summary = summarize_results(resFile)
    % Provide details from results file
    summary = struct;
    data = load(resFile);
    results = data.results;
    % Total number of samples examined
    summary.N = length(results);
    % Verified (unsat)
    summary.unsat = sum(results(:,1) == 1);
    % Falsified (sat) from input set
    summary.sat = sum(results(:,1) == 0);
    % Unknown (not using exact)
    summary.unknown = sum(results(:,1) == 2);
    % Misclasified (sat) of original image
    summary.missclassified = sum(results(:,1) == -1);
    % Also possible -2 (error), which may be out of memory (most common here)
    summary.avgTime = sum(results(:,2))/summary.N;
end