%% Get all results, one dataset at a time

results = dir('results');

% Variable in the study
% datasets = ["adrenal", "fracture", "nodule", "organ", "synapse", "vessel"];
datasets = ["fracture", "nodule", "organ", "synapse"]; 
attackType = ["bright", "dark"];
maxpixels = ["50", "100", "500", "1000"];
epsilon = ["2", "5", "10", "25"]; % epsilon/255

% Initialize directories
mkdir('plots');
for i=datasets
    mkdir("plots/"+i);
end

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
 

% Time for plots (maxpixels variable)
for ds = datasets
    for adv = attackType
        for ep = epsilon
            
            % Initialize vars to plot
            sat   = [];
            unsat = [];
            unk   = [];
            miss  = [];
            avgVT = [];
            avgRT = [];
            
            % Get data
            for mp = maxpixels
                resFile = "results/verification_" + ds + "_" + adv + "_" + ep + "_" + mp +".mat";
                res = summarize_results(resFile);
                sat = [sat, res.sat];
                unsat = [unsat, res.unsat];
                unk = [unk, res.unknown];
                miss = [miss, res.misclassified];
                avgVT = [avgVT, res.avgTime]; % average computation time to verify one instance
                avgRT = [avgRT, res.avgRT]; % average computation time to compute the reachable sets (miss and sat do not count here)
            end
            counts = [unsat; sat; unk; miss];
            
            
            % Create figure (results)
            f = figure('visible','off');
            bar(1:4, counts','stacked') % plot verification results
            grid;
            % set values for x-axis
            xticks([1 2 3 4]);
            xticklabels(maxpixels);
            % axis labels
            xlabel('Max pixels perturbed');
            ylabel('# instances')
            % legend
            leg = {"unsat", "sat", "unknown", "missclass"};
            legend(leg, 'Location', 'best');
            % save figure
            saveas(f, "plots/" + ds + "/verification_" + adv + "_" + ep + ".png");
            
            % create figure (time)
            f = figure('visible','off');
            % plot computation time results
            plot(1:4, avgVT, 'r--o');
            hold on;
            plot(1:4, avgRT, 'b--v');
            % set values for x-axis
            xticks([1 2 3 4]);
            xticklabels(maxpixels)
            % axis labels
            xlabel('Max pixels perturbed');
            ylabel("Average Time (s)")
            % legend
            leg = {"all", "unsat & unknown"};
            legend(leg, 'Location', 'best');
            % save figure
            saveas(f, "plots/" + ds + "/avgTime_" + adv + "_" + ep + ".png");

        end
    end
end

% Time for plots (epsilon variable)
for ds = datasets
    for adv = attackType
        for mp = maxpixels
            
            % Initialize vars to plot
            sat   = [];
            unsat = [];
            unk   = [];
            miss  = [];
            avgVT = [];
            avgRT = [];
            
            % Get data
            for ep = epsilon
                resFile = "results/verification_" + ds + "_" + adv + "_" + ep + "_" + mp +".mat";
                res = summarize_results(resFile);
                sat = [sat, res.sat];
                unsat = [unsat, res.unsat];
                unk = [unk, res.unknown];
                miss = [miss, res.misclassified];
                avgVT = [avgVT, res.avgTime]; % average computation time to verify one instance
                avgRT = [avgRT, res.avgRT]; % average computation time to compute the reachable sets (miss and sat do not count here)
            end
            counts = [unsat; sat; unk; miss];
            
            % Create figure (results)
            f = figure('visible','off');
            bar(1:3, counts','stacked') % plot results
            grid
            % set values for x-axis
            xticks([1 2 3]);
            xticklabels(epsilon)
            % axis labels
            xlabel('epsilon (\epsilon)');
            ylabel('# instances')
            % legend
            leg = {"unsat", "sat", "unknown", "missclass"};
            legend(leg, 'Location', 'best');
            % save figure
            saveas(f, "plots/" + ds + "/verification_" + adv + "_" + mp + ".png");
            
            % create figure (time)
            f = figure('visible','off');
            % plot timing results
            plot(1:3, avgVT, 'r--o');
            hold on;
            plot(1:3, avgRT, 'b--v');
            % set values for x-axis
            xticks([1 2 3]);
            xticklabels(epsilon)
            % axis labels
            xlabel('epsilon (\epsilon)');
            ylabel("Average Time (s)")
            % legend
            leg = {"all", "unsat & unknown"};
            legend(leg, 'Location', 'best');
            % save figure
            saveas(f, "plots/" + ds + "/avgTime_" + adv + "_" + mp + ".png");

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
    
    % Find instances where reachability is needed
    x = [find(results(:,1)==2); find(results(:,1)==1)];
    nx = length(x); % how many instances reachability is computed for

    % Average time only for unknown and unsat properties (reachsets computed)
    summary.avgRT = sum(results(x,2))/nx;
end