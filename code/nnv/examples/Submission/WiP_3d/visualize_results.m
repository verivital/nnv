%% Get all results, one dataset at a time

results = dir('results');

% Variable in the study
% datasets = ["adrenal", "fracture", "nodule", "organ", "synapse", "vessel"];
datasets = ["fracture", "nodule", "organ", "synapse"]; 
attackType = ["bright", "dark"];
maxpixels = ["50", "100", "500", "1000"];
epsilon = ["2", "4", "10"]; % epsilon/255

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
 

% Time for plots
for ds = datasets
    for adv = attackType
        for ep = epsilon
            % Initialize vars to plot
            sat   = [];
            unsat = [];
            unk   = [];
            miss  = [];
            avgTime = [];
            % Get data
            for mp = maxpixels
                resFile = "verification_" + ds + "_" + adv + "_" + ep + "_" + mp +".mat";
                res = summarize_results(resFile);
                sat = [sat, res.sat];
                unsat = [unsat, res.unsat];
                unk = [unk, res.unknown];
                miss = [miss, res.misclassified];
                avgTime = [avgTime, res.avgTime];
            end
            counts = [unsat; sat; unk; miss];
            leg = {"unsat", "sat", "unknown", "missclass"};
            % Create figure
            f = figure;
            bar(1:4, counts','stacked') % change ticks and label later on
            grid
            legend(leg, 'Location', 'best');
            saveas(f, "plots/verification_" + ds + "_" + adv + "_" + ep+".png");
        end
    end
end


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
    summary.misclassified = sum(results(:,1) == -1);
    % Also possible -2 (error), which may be out of memory (most common here)
    summary.avgTime = sum(results(:,2))/summary.N;
end