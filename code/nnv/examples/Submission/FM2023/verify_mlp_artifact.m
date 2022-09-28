%% MAKE SURE YOU ARE INSIDE THE FM2023/ FOLDER

format long

%% MLP MODEL ID: INPUT AN ID FROM 0 TO 4 TO CHOOSE THE MODEL (LINE 10)
% 0 - MLP0
% 1 - MLP1
% 2 - MLP2
% 3 - MLP3
% 4 - MLP4
current_model = 0;

%% STAR REACHABILITY APPROACH: LINES 15,16 CONTAIN THE REACHABILILTY METHOD OPTIONS
% 'exact-star' - exact reachability
% 'approx-star' - over-approximate reachability
% approach = 'exact-star';
approach = 'approx-star';
    
if strcmp(approach, 'exact-star')
    approachName = 'exact';
else
    approachName = 'approx';
end

%% MLP DATA PATHS
dataFolder = strcat("mlp_comparison/Marabou/digit_test_set.csv");

data = load(dataFolder);
    
X = data(:, 2:size(data,2));

if current_model == 0 % mlp_init
    bnn = load('mlp_comparison/nn_models/mlp0.mat').bnn;
    model_name = "mlp0";
elseif current_model == 1 % mlp8
    bnn = load('mlp_comparison/nn_models/mlp1.mat').bnn;
    path = "mlp_comparison/mlp1_4_examples/mlp1_ids.mat";
    model_name = "mlp1";
elseif current_model == 2 % mlp9
    bnn = load('mlp_comparison/nn_models/mlp2.mat').bnn;
    path = "mlp_comparison/mlp1_4_examples/mlp2_ids.mat";
    model_name = "mlp2";
elseif current_model == 3 % mlp10
    bnn = load('mlp_comparison/nn_models/mlp3.mat').bnn;
    path = "mlp_comparison/mlp1_4_examples/mlp3_ids.mat";
    model_name = "mlp3";
elseif current_model == 4
    bnn = load('mlp_comparison/nn_models/mlp4.mat').bnn;
    path = "mlp_comparison/mlp1_4_examples/mlp4_ids.mat";
    model_name = "mlp4";
end


mpl_results_folder = "mlp_comparison/Marabou/mlp_results.csv";

mlp_results_data = readtable(mpl_results_folder);

if current_model == 0
    sample_data = zeros(size(mlp_results_data, 1), 3);
    
    for i=1:size(mlp_results_data, 1)    
        if i > 300 && i < 401
            sample_data(i, 1) = str2double(extractBetween(mlp_results_data.Var1{i, 1}, strfind(mlp_results_data.Var1{i, 1}, '.') - 2, strfind(mlp_results_data.Var1{i, 1}, '.') + 1));
            sample_data(i, 2) = str2double(extractBetween(mlp_results_data.Var1{i, 1}, length(mlp_results_data.Var1{i, 1}) - 3, length(mlp_results_data.Var1{i, 1})));
        elseif i > 50 && i < 101
            sample_data(i, 1) = str2double(extractBetween(mlp_results_data.Var1{i, 1}, strfind(mlp_results_data.Var1{i, 1}, '.') - 1, strfind(mlp_results_data.Var1{i, 1}, '.') + 2));
            sample_data(i, 2) = str2double(extractBetween(mlp_results_data.Var1{i, 1}, length(mlp_results_data.Var1{i, 1}) - 3, length(mlp_results_data.Var1{i, 1})));
        else
            sample_data(i, 1) = str2double(extractBetween(mlp_results_data.Var1{i, 1}, strfind(mlp_results_data.Var1{i, 1}, '.') - 1, strfind(mlp_results_data.Var1{i, 1}, '.') + 1));
            sample_data(i, 2) = str2double(extractBetween(mlp_results_data.Var1{i, 1}, length(mlp_results_data.Var1{i, 1}) - 3, length(mlp_results_data.Var1{i, 1})));
        end

        if mlp_results_data.Var2{i, 1} == "unsat"
            sample_data(i, 3) = 1;
        elseif mlp_results_data.Var2{i, 1} == "sat"
            sample_data(i, 3) = 0;
        else
            sample_data(i, 3) = -1; 
        end
    end
else
    sample_data = zeros(size(500, 1), 3);

    ids = load(path).ids;

    for i=1:length(deltas)
        for j=1:50
            sample_data((i - 1) * 50 + j, 1) = deltas(i);
            sample_data((i - 1) * 50 + j, 2) = ids((i - 1) * 50 + j);
            sample_data((i - 1) * 50 + j, 3) = -1;
        end
    end
end


verdict_star = zeros(size(sample_data, 1), 1);

corrected_samples_num = 0;
setback_samples_num = 0;

reachT = {};
signReachT = {};

deltas = [0.1, 0.15, 0.2, 0.3, 0.5, 1, 3, 5, 10, 15];

starUnsatMap = containers.Map(deltas, zeros(1, size(deltas, 2)));
starSatMap = containers.Map(deltas, zeros(1, size(deltas, 2)));
starUnkMap = containers.Map(deltas, zeros(1, size(deltas, 2)));

reluplexUnsatMap = containers.Map(deltas, zeros(1, size(deltas, 2)));
reluplexSatMap = containers.Map(deltas, zeros(1, size(deltas, 2)));
reluplexUnkMap = containers.Map(deltas, zeros(1, size(deltas, 2)));

for i=1:size(sample_data, 1)    
    class = data(sample_data(i, 2), 1);
    
    %% extra eval mlp
    img = X(sample_data(i, 2), :)';
    
    lb = img - sample_data(i, 1);
    ub = img + sample_data(i, 1);

    lb(lb < 0) = 0;
    lb(lb > 255) = 255;
    ub(ub < 0) = 0;
    ub (ub > 255) = 255;

    S = Star(lb, ub);
    
    %% Verify
    [~, reachT{i}, signReachT{i}] = bnn.reach(S, approach, 1);
    
    [val, pos] = max(bnn.outputSet.V(:,1));
    
    eval_vector = bnn.evaluate(img);
    [res, evRes] = max(eval_vector);
    
    if class == evRes - 1
        if evRes == pos
            if sum(bnn.outputSet.V(:,1) == val) == 1
                verdict_star(i) = 1;

                starUnsatMap(sample_data(i, 1)) = starUnsatMap(sample_data(i, 1)) + 1;

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
            end
          else
            if strcmp(approach, 'exact-star')
                verdict_star(i) = 0;

                starSatMap(sample_data(i, 1)) = starSatMap(sample_data(i, 1)) + 1;

                if sample_data(i, 3) == 1
                    reluplexUnsatMap(sample_data(i, 1)) = reluplexUnsatMap(sample_data(i, 1)) + 1;
                else
                    if sample_data(i, 3) == 0 
                        reluplexSatMap(sample_data(i, 1)) = reluplexSatMap(sample_data(i, 1)) + 1;
                    else
                        reluplexUnkMap(sample_data(i, 1)) = reluplexUnkMap(sample_data(i, 1)) + 1;
                    end
                end
            elseif strcmp(approach, 'approx-star')
                [flag, counterEx] = findCounterEx(S, bnn, 50, evRes); % find a counterexample if SAT is returned by
                                                                      % over-approximate
                                                                      % analysis

                if flag % counterexample was found => SAT
                    verdict_star(i) = 0;

                    starSatMap(sample_data(i, 1)) = starSatMap(sample_data(i, 1)) + 1;

                    if sample_data(i, 3) == 1
                        reluplexUnsatMap(sample_data(i, 1)) = reluplexUnsatMap(sample_data(i, 1)) + 1;
                    else
                        if sample_data(i, 3) == 0 
                            reluplexSatMap(sample_data(i, 1)) = reluplexSatMap(sample_data(i, 1)) + 1;
                        else
                            reluplexUnkMap(sample_data(i, 1)) = reluplexUnkMap(sample_data(i, 1)) + 1;
                        end
                    end
                else % counterexample wasn't found => UNK
                    verdict_star(i) = -1;
                    starUnkMap(sample_data(i, 1)) = starUnkMap(sample_data(i, 1)) + 1;

                    if sample_data(i, 3) == 1
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
    end
end

sample_data(:, 4) = verdict_star;

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

save(strcat('mlp_comparison/results/', model_name, '_', approachName,'_result.mat'), 'sample_data');
save(strcat('mlp_comparison/results/', model_name, '_', approachName,'_reachT.mat'), 'reachT');
save(strcat('mlp_comparison/results/', model_name, '_', approachName,'_signReachT.mat'), 'signReachT');

%% Functions:

function [flag, counterEx] = findCounterEx(In, bnn, samplesNum, groundLabel)
    flag = false;
    counterEx = [];

    Out = bnn.outputSet;
    cL = size(In.V, 1);
    l = size(In.V(:,2:end), 2);

    flag = false;
    
    inputs = zeros(cL, samplesNum);
    
    for i = 1:samplesNum
    
    %% init mlp
    %S = S(i);
        predicateVecs = zeros(cL, l);
        inputs(:, i) = In.V(:, 1);
        
        for j = 1:l
            lb = Out.predicate_lb(j);
            ub = Out.predicate_ub(j);
            
            predicateVecs(:, j) = lb + (ub - lb) .* rand(cL,1);
            inputs(:, i) = inputs(:, i) + In.V(:, j + 1) .* predicateVecs(:, j);
        end
                
        [res, pos] = max(bnn.evaluate(inputs(:, i)));
        
        if pos ~= groundLabel
            flag = true;
            counterEx = inputs(:, i);
            return;
        end
    end
end

function S = constructInput(data, disturbances)
    S = [];

    for i=1:size(data, 2)
        lb = data(:, i) - disturbances(i, 1);
        ub = data(:, i) + disturbances(i, 1);

        lb(lb < 0) = 0;
        lb(lb > 255) = 255;
        ub(ub < 0) = 0;
        ub (ub > 255) = 255;

        S = [S Star(lb, ub)];
    end
end

function v = constructVec(data)
    v = zeros(size(data,3),1);
    
    for i = 1:length(data)
        v(i) = data(1,1,i);
    end
end