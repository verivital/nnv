%% MAKE SURE YOU ARE INSIDE THE CAV_2022/ FOLDER

format long

bnn = load('xnor_comparison/nn_models/xnor0.mat').bnn;

X = load('xnor_comparison/Marabou/fashion_test_set.csv');
    
xnor_results_folder = "xnor_comparison/Marabou/xnor_results.csv";

xnor_results_data = readtable(xnor_results_folder);

sample_data = zeros(size(xnor_results_data, 1), 3);

xnor_class_ids = [];

for i=1:size(xnor_results_data, 1)
    sample_data(i, 1) = xnor_results_data.Var6(i);
    a = strfind(xnor_results_data.Var9{i, 1}, ',');
    sample_data(i, 2) = str2double(extractBetween(xnor_results_data.Var9{i, 1}, 1, a(1) - 1));

    if strfind(xnor_results_data.Var9{i, 1}, 'unsat') > 0
        sample_data(i, 3) = 1;
    elseif strfind(xnor_results_data.Var9{i, 1}, 'sat') > 0
        sample_data(i, 3) = 0;
    else
        sample_data(i, 3) = -1; 
    end

    [res, pos] = max(bnn.evaluate(reshape(X(sample_data(i,2), :), [28,28,1])'));
end


%% RESULTS CONTAINERS
verdict_star = zeros(size(sample_data, 1), 1);

corrected_samples_num = 0;
setback_samples_num = 0;

reachT = {};
signReachT = {};

deltas = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3];

starUnsatMap = containers.Map(deltas, zeros(1, size(deltas, 2)));
starSatMap = containers.Map(deltas, zeros(1, size(deltas, 2)));
starUnkMap = containers.Map(deltas, zeros(1, size(deltas, 2)));


reluplexUnsatMap = containers.Map(deltas, zeros(1, size(deltas, 2)));
reluplexSatMap = containers.Map(deltas, zeros(1, size(deltas, 2)));
reluplexUnkMap = containers.Map(deltas, zeros(1, size(deltas, 2)));

%% VERIFICATION
for i=1:size(sample_data, 1)
    img = reshape(X(sample_data(i, 2), :), [28, 28, 1])';
    evRes = bnn.evaluate(img);
    
    % get the ground-truth class (for mlp0 all the given examples by
    % Artifact are claimed to be classified correctly
    [~, class] = max(evRes);
    
    % construct a Star
    lb = img - sample_data(i, 1);
    ub = img + sample_data(i, 1);

    lb(lb < 0) = 0;
    lb(lb > 255) = 255;
    ub(ub < 0) = 0;
    ub (ub > 255) = 255;

    S = ImageStar(lb, ub);
        
    %% Verify
    [~, reachT{i}, signReachT{i}] = bnn.reach(S, 'approx-star', 1);
    
    [val, pos] = max(bnn.outputSet.V(:,1));
    
    if class == pos
        % Star classifies as UNSAT
        if sum(bnn.outputSet.V(:,1) == val) == 1
            verdict_star(i) = 1;
            
            starUnsatMap(sample_data(i, 1)) = starUnsatMap(sample_data(i, 1)) + 1;

            % Marabou classifies as SAT or UNKNOWN
            if sample_data(i, 3) < 1 
                corrected_samples_num = corrected_samples_num + 1;

                if sample_data(i, 3) == 0 
                    reluplexSatMap(sample_data(i, 1)) = reluplexSatMap(sample_data(i, 1)) + 1;
                else
                    reluplexUnkMap(sample_data(i, 1)) = reluplexUnkMap(sample_data(i, 1)) + 1;
                end
            else
                reluplexUnsatMap(sample_data(i, 1)) = reluplexUnsatMap(sample_data(i, 1)) + 1;
            end
        % Star classifies as SAT
        else
            verdict_star(i) = -1;
            
            starUnkMap(sample_data(i, 1)) = starUnkMap(sample_data(i, 1)) + 1;

            % Marabou classifies as UNSAT
            if sample_data(i, 3) == 1
                setback_samples_num = setback_samples_num + 1;

                reluplexUnsatMap(sample_data(i, 1)) = reluplexUnsatMap(sample_data(i, 1)) + 1;
            else
                if sample_data(i, 3) == 0 
                    reluplexSatMap(sample_data(i, 1)) = reluplexSatMap(sample_data(i, 1)) + 1;
                else
                    reluplexUnkMap(sample_data(i, 1)) = reluplexUnkMap(sample_data(i, 1)) + 1;
                end
            end
        end
    % Star can't decide between several classes (relatively equal probabilities -> SAT)    
    else
        [flag, counterEx] = findCounterEx(S, bnn, 200, class);

        if flag
            verdict_star(i) = 0;

            starSatMap(sample_data(i, 1)) = starSatMap(sample_data(i, 1)) + 1;

            if sample_data(i, 3) == 1
                setback_samples_num = setback_samples_num + 1;
                
                reluplexUnsatMap(sample_data(i, 1)) = reluplexUnsatMap(sample_data(i, 1)) + 1;
            else
                if sample_data(i, 3) == 0 
                    reluplexSatMap(sample_data(i, 1)) = reluplexSatMap(sample_data(i, 1)) + 1;
                else
                    reluplexUnkMap(sample_data(i, 1)) = reluplexUnkMap(sample_data(i, 1)) + 1;
                end
            end
        else
            verdict_star(i) = -1;
            starUnkMap(sample_data(i, 1)) = starUnkMap(sample_data(i, 1)) + 1;

            if sample_data(i, 3) == 1
                setback_samples_num = setback_samples_num + 1;
                
                reluplexUnsatMap(sample_data(i, 1)) = reluplexUnsatMap(sample_data(i, 1)) + 1;
            else
                if sample_data(i, 3) == 0 
                    reluplexSatMap(sample_data(i, 1)) = reluplexSatMap(sample_data(i, 1)) + 1;
                else
                    reluplexUnkMap(sample_data(i, 1)) = reluplexUnkMap(sample_data(i, 1)) + 1;
                end
            end
        end
    end
end

sample_data(:, 4) = verdict_star;

time_data = zeros(size(deltas,1), 3);

%% Print results:
fprintf("\nThe number of UNSAT samples (Reluplex-Star):\t\t\t\t%d\t---\t%d\n", length(find(sample_data(:, 3) == 1)), length(find(sample_data(:, 4) == 1)));

fprintf("The number of SAT samples (Reluplex-Star):\t\t\t\t%d\t---\t%d\n", length(find(sample_data(:, 3) == 0)), length(find(sample_data(:, 4) == 0)));

fprintf("The number of unsolved samples (Reluplex-Star):\t\t\t\t%d\t---\t%d\n", length(find(sample_data(:, 3) == -1)), length(find(sample_data(:, 4) == -1)));

for i=1:size(deltas, 2)
    fprintf("\nThe number of UNSAT samples for delta %.5f (Reluplex-Star):\t\t%d\t---\t%d\n", deltas(i), reluplexUnsatMap(deltas(i)), starUnsatMap(deltas(i)));
    fprintf("The number of SAT samples for delta %.5f (Reluplex-Star):\t\t%d\t---\t%d\n", deltas(i), reluplexSatMap(deltas(i)), starSatMap(deltas(i)));
    fprintf("The number of UNK samples for delta %.5f (Reluplex-Star):\t\t%d\t---\t%d\n", deltas(i), reluplexUnkMap(deltas(i)), starUnkMap(deltas(i)));
    
    unsat_ids = find(sample_data(:,1) == deltas(i) & sample_data(:,4) == 1);
    sat_ids = find(sample_data(:,1) == deltas(i) & sample_data(:,4) == 0);
    
    fprintf("Avg. UNSAT reach time for delta %.5f (Reluplex-Star):\t\t%d\n", deltas(i), sum([reachT{unsat_ids}]) / length([reachT{unsat_ids}]));
    fprintf("Avg. SAT reach time for delta %.5f (Reluplex-Star):\t\t\t%d\n", deltas(i), sum([reachT{sat_ids}]) / length([reachT{sat_ids}]));
end

fprintf("\nThe general reach time:\t\t%d\n", sum([reachT{:}]));
fprintf("\nThe general sign time:\t\t%d\n", sum([signReachT{:}]));

fprintf("Avg reach time:\t\t\t%d\n", mean([reachT{:}]));
fprintf("Avg SignLayer reach time:\t%d\n", mean([signReachT{:}]));

save(strcat('xnor_comparison/results/xnor0_approx_result.mat'), 'sample_data');
save(strcat('xnor_comparison/results/xnor0_approx_reachT.mat'), 'reachT');
save(strcat('xnor_comparison/results/xnor0_approx_signReachT.mat'), 'signReachT');

function [flag, counterEx] = findCounterEx(I, bnn, samplesNum, groundLabel)
    counterEx = [];
    Out = bnn.outputSet;
    In = I.toStar();
    cL = size(In.V, 1);
    l = size(In.V(:,2:end), 2);

    flag = false;
    
    inputs = zeros(cL, samplesNum);
    
    for i = 1:samplesNum
        predicateVecs = zeros(cL, l);
        inputs(:, i) = In.V(:, 1);
        
        for j = 1:l
            lb = Out.predicate_lb(j);
            ub = Out.predicate_ub(j);
            
            predicateVecs(:, j) = lb + (ub - lb) .* rand(cL,1);
            inputs(:, i) = inputs(:, i) + In.V(:, j + 1) .* predicateVecs(:, j);
        end
        
        currentInput = reshape(inputs(:, i), 28, 28);
        
        [res, pos] = max(bnn.evaluate(currentInput));
        
        if pos ~= groundLabel
            flag = true;
            counterEx = inputs(:, i);
            return;
        end
    end
end