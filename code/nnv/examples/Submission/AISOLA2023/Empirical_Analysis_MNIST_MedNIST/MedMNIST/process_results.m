%% Process results
results_folder = dir('results');

% N = 45;
regs = ["dropout", "jacobian", "l2"];
Nregs = length(regs);
inits = ["glorot", "he", "narrow"]; 
Ninits = length(inits);
abdomen = 1:50;
breast  = 51:100;
chest   = 101:150;
cxr     = 151:200;
hand    = 201:250;
head    = 251:300;
classes = [abdomen; breast; chest; cxr; hand; head];
numClasses = 6;

% first load all results
N = length(results_folder)-2;
allRes = cell(N,2);

for i=3:N+2
    ind_res = load(['results', filesep, results_folder(i).name]);
    allRes{i-2,1} = ind_res.res;
    allRes{i-2,2} = results_folder(i).name;
end

% results are ordered from 
% - dropout -> jacobian -> l2
% - within each, we have: glorot -> he -> narrow normal
% - within each, 0 -> 1 -> 2 -> 3-> 4

%% Analysis per class

% res indexes per class
% abdomen = 1:50;
% breast  = 51:100;
% chest   = 101:150;
% cxr     = 151:200;
% hand    = 201:250;
% head    = 251:300;
% 
% classes = [abdomen; breast; chest; cxr; hand; head];
% numClasses = 6;
classRes = zeros(numClasses,4);

% Process results per class
for c = 1:numClasses
    idxs = classes(c,:);
    Rob = 0; Unk = 0; Norob = 0; Avgtime = 0;
    for i=1:N
        [rob, unk, norob, time] = process_model_res(allRes{i,1}, idxs);
        Rob = Rob + rob; Unk = Unk + unk; 
        Norob = Norob + norob; Avgtime = Avgtime + time;
    end
    Avgtime = Avgtime/N;
    classRes(c,:) = [Rob, Unk, Norob, Avgtime];
end

%% Analysis per regularizer

% regs = ["dropout", "jacobian", "l2"];
% Nregs = length(regs);
% regRes = zeros(Nregs,4); % dropout, jacobian, l2
% idxs = 1:300; % all indexes for the rest of the results
% 
% for r = 1:Nregs
%     Rob = 0; Unk = 0; Norob = 0; Avgtime = 0;
%     count = 0;
%     for i=1:N
%         if contains(allRes{i,2}, regs(r))
%             [rob, unk, norob, time] = process_model_res(allRes{i,1}, idxs);
%             Rob = Rob + rob; Unk = Unk + unk; 
%             Norob = Norob + norob; Avgtime = Avgtime + time;
%             count = count + 1;
%         end
%     end
%     Avgtime = Avgtime/count;
%     regRes(r,:) = [Rob, Unk, Norob, Avgtime];
% end

%% Analysis per initialization scheme

% inits = ["glorot", "he", "narrow"]; 
% Ninits = length(inits);
% initRes = zeros(Ninits,4);
% idxs = 1:300; % all indexes for the rest of the results
% 
% for r = 1:Ninits
%     Rob = 0; Unk = 0; Norob = 0; Avgtime = 0;
%     count = 0;
%     for i=1:N
%         if contains(allRes{i,2}, inits(r))
%             [rob, unk, norob, time] = process_model_res(allRes{i,1}, idxs);
%             Rob = Rob + rob; Unk = Unk + unk; 
%             Norob = Norob + norob; Avgtime = Avgtime + time;
%             count = count + 1;
%         end
%     end
%     Avgtime = Avgtime/count;
%     initRes(r,:) = [Rob, Unk, Norob, Avgtime];
% end

%% Analysis per seed (0,1,2,3,4)

seeds = ["0", "1", "2", "3", "4"]; 
Nseeds = length(seeds);
seedRes = zeros(Nseeds,4);
idxs = 1:300; % all indexes for the rest of the results

% Compute average across all models
for i=1:N
    loc_seed = mod(i,5);
    if loc_seed == 0
        loc_seed = 5;
    end
    [rob, unk, norob, time] = process_model_res(allRes{i,1}, idxs);
    seedRes(loc_seed, 1) = seedRes(loc_seed, 1)  + rob/2700; % 9 (3 init x 3 reg) * 300 (images) = 2700
    seedRes(loc_seed, 2) = seedRes(loc_seed, 2)  + unk/2700;
    seedRes(loc_seed, 3) = seedRes(loc_seed, 3)  + norob/2700;
    seedRes(loc_seed, 4) = seedRes(loc_seed, 4)  + time/9;
end

% Individually
% Seed 0
model_idxs_0 = [1,6,11,16,21,26,31,36,41]';
stats_seed0 = process_model_classes(allRes, classes, model_idxs_0);

% seed 1
model_idxs_1 = model_idxs_0 + 1;
stats_seed1 = process_model_classes(allRes, classes, model_idxs_1);

% seed 2
model_idxs_2 = model_idxs_1 + 1;
stats_seed2 = process_model_classes(allRes, classes, model_idxs_2);

% seed 3
model_idxs_3 = model_idxs_2 + 1;
stats_seed3 = process_model_classes(allRes, classes, model_idxs_3);

% seed 4
model_idxs_4 = model_idxs_3 + 1;
stats_seed4 = process_model_classes(allRes, classes, model_idxs_4);


%% Analysis per reg, per init (3 x 3)

% Order:
% dropout (glorot -> he -> narrow) -> jacobian (glorot -> he -> narrow) -> L2 (glorot -> he -> narrow)

% regInitRes = zeros(Ninits*Nregs,4);
% idxs = 1:300; % all indexes for the rest of the results
% 
% for i=1:N
%     combo = ceil(i/5);
%     [rob, unk, norob, time] = process_model_res(allRes{i,1}, idxs);
%     regInitRes(combo,1) = regInitRes(combo,1) + rob; % robust
%     regInitRes(combo,2) = regInitRes(combo,2) + unk; % unknown
%     regInitRes(combo,3) = regInitRes(combo,3) + norob; % not robust
%     regInitRes(combo,4) = regInitRes(combo,4) + time; % computation time
% end
% regInitRes(:,4) = regInitRes(:,4)/15; % 5 = number of models per combo


%% Analysis per reg per class

% Order:
% dropout (1,2,3,4,5,6) -> jacobian (1,2,3,4,5,6) -> L2 (1,2,3,4,5,6)
% 
% regClassRes = zeros(numClasses*Nregs,4);
% 
% % Process results per class
% for c = 1:numClasses
%     idxs = classes(c,:);
%     for i=1:N
%         [rob, unk, norob, time] = process_model_res(allRes{i,1}, idxs);
%         if contains(allRes{i,2}, "dropout")
%             regClassRes(c,1) = regClassRes(c,1) + rob;
%             regClassRes(c,2) = regClassRes(c,2) + unk;
%             regClassRes(c,3) = regClassRes(c,3) + norob;
%             regClassRes(c,4) = regClassRes(c,4) + time;
%         elseif contains(allRes{i,2}, "jacobian")
%             regClassRes(c+6,1) = regClassRes(c+6,1) + rob;
%             regClassRes(c+6,2) = regClassRes(c+6,2) + unk;
%             regClassRes(c+6,3) = regClassRes(c+6,3) + norob;
%             regClassRes(c+6,4) = regClassRes(c+6,4) + time;
%         else % L2
%             regClassRes(c+12,1) = regClassRes(c+12,1) + rob;
%             regClassRes(c+12,2) = regClassRes(c+12,2) + unk;
%             regClassRes(c+12,3) = regClassRes(c+12,3) + norob;
%             regClassRes(c+12,4) = regClassRes(c+12,4) + time;
%         end
%     end
% end
% regClassRes(:,4) = regClassRes(:,4)/15; % (5 models per init, so 15 models total)


%% Analysis per reg per init per class

regInitClassRes = zeros(numClasses*Nregs*Ninits,4);
% dropout (1 to 18)
% jacobian (19 to 36)
% l2 (37 to 54)

for i=1:N
    combo = floor((i-1)/5);
    for c = 1:numClasses
        idxs = classes(c,:);
        [rob, unk, norob, time] = process_model_res(allRes{i,1}, idxs);
        regInitClassRes(combo*6+c,1) = regInitClassRes(combo*6+c,1) + rob; % robust
        regInitClassRes(combo*6+c,2) = regInitClassRes(combo*6+c,2) + unk; % unknown
        regInitClassRes(combo*6+c,3) = regInitClassRes(combo*6+c,3) + norob; % not robust
        regInitClassRes(combo*6+c,4) = regInitClassRes(combo*6+c,4) + time; % computation time
    end
end
regInitClassRes(:,4) = regInitClassRes(:,4)/5; % 5 = number of models per combo

classNames = ["AbdomenCT", "BreastMRI", "ChestCT", "CXR", "Hand", "HeadCT"];


%% All stats for each model
all_idxs = 1:45;
all_stats = process_model_classes(allRes, classes, all_idxs');
summary_res = sum(all_stats,2);
summary_res = squeeze(summary_res);
summary_res(:,4) = summary_res(:,4)/numClasses;
% worst models 8 (203/300), 10 (196/300), 38 (196/300) and 40 (188/300)
% 8,10 -> dropout_he
% 38,40 -> l2_he
% slowest models are 6 (43.8), 9 (42.1), 10 (43.9), 39 (42.45) and 40 (44.4)
% 6,9,10 -> dropout_he
% 39,40 -> l2_he
% worst model is **40**
%
% best models are -> 11, 15 (dropout_narrow), 26, 28, 29, and 30 (jacobian_narrow)
% fastest model -> 14 (dropout_narrow) with 5 seconds, significantly faster than others

%% Visualize results


%% Compare regularization vs classes

dropout_class = regInitClassRes(1:6,:) + regInitClassRes(7:12,:) + regInitClassRes(13:18,:);
dropout_class(:,4) = dropout_class(:,4)/3; % compute average time
dropout_class(:,1:3) = dropout_class(:,1:3)/750; % total number of instances (plot as percentage)
jacobian_class = regInitClassRes(19:24,:) + regInitClassRes(25:30,:) + regInitClassRes(31:36,:);
jacobian_class(:,4) = jacobian_class(:,4)/3;
jacobian_class(:,1:3) = jacobian_class(:,1:3)/750; % total number of instances (plot as percentage)
l2_class= regInitClassRes(37:42,:) + regInitClassRes(43:48,:) + regInitClassRes(49:54,:);
l2_class(:,4) = l2_class(:,4)/3;
l2_class(:,1:3) = l2_class(:,1:3)/750; % total number of instances (plot as percentage)

% Create figure
figure;
grid;hold on;
plot(1:6, dropout_class(:,1),'r-d');
plot(1:6, jacobian_class(:,1),'b-o');
plot(1:6, l2_class(:,1),'k--');
set(gca, 'xtick', 1:6)
set(gca, 'xticklabel', classNames);
% set(gca, "YTick", 0.9:0.01:1);
% ylim([0.925, 1.005])
ylabel("Robust %");
legend('dropout','jacobian', 'L2', 'Location','best');
exportgraphics(gca, "plots/regRes_vs_class.pdf",'ContentType','vector');

% Create figure
figure;
grid;hold on;
plot(1:6, dropout_class(:,4),'r-d');
plot(1:6, jacobian_class(:,4),'b-o');
plot(1:6, l2_class(:,4),'k--');
set(gca, 'xtick', 1:6)
set(gca, 'xticklabel', classNames);
ylabel("Time (s)")
legend('dropout','jacobian', 'L2', 'Location','best');
exportgraphics(gca, "plots/regTime_vs_class.pdf",'ContentType','vector');


%% Compare initializers vs classes 

glorot_class = regInitClassRes(1:6,:) + regInitClassRes(19:24,:) + regInitClassRes(37:42,:);
glorot_class(:,4) = glorot_class(:,4)/3; % compute average time
glorot_class(:,1:3) = glorot_class(:,1:3)/750; % total number of instanes (plot as percentage)
he_class = regInitClassRes(7:12,:) + regInitClassRes(25:30,:) + regInitClassRes(43:48,:) ;
he_class(:,4) = he_class(:,4)/3;
he_class(:,1:3) = he_class(:,1:3)/750; % total number of instanes (plot as percentage)
narrow_class = regInitClassRes(13:18,:) + regInitClassRes(31:36,:) + regInitClassRes(49:54,:);
narrow_class(:,4) = narrow_class(:,4)/3;
narrow_class(:,1:3) = narrow_class(:,1:3)/750; % total number of instanes (plot as percentage)

% Create figure
figure;
grid;hold on;
plot(1:6, glorot_class(:,1),'r-d');
plot(1:6, he_class(:,1),'b-o');
plot(1:6, narrow_class(:,1),'k--');
set(gca, 'xtick', 1:6)
set(gca, 'xticklabel', classNames);
% set(gca, "YTick", 0.9:0.01:1);
% ylim([0.925, 1.005])
ylabel("Robust %");
legend('glorot','he', 'narrow-normal', 'Location','best');
exportgraphics(gca, "plots/initRes_vs_class.pdf",'ContentType','vector');

% Create figure
figure;
grid;hold on;
plot(1:6, glorot_class(:,4),'r-d');
plot(1:6, he_class(:,4),'b-o');
plot(1:6, narrow_class(:,4),'k--');
set(gca, 'xtick', 1:6)
set(gca, 'xticklabel', classNames);
ylabel("Time (s)")
legend('glorot','he', 'narrow-normal', 'Location','best');
exportgraphics(gca, "plots/initTime_vs_class.pdf",'ContentType','vector');

%% Compare combinations vs classes 

% Create figure
figure;
grid;hold on;
plot(1:6, regInitClassRes(1:6,1)/250,'--d', 'Color', "#A2142F");
plot(1:6, regInitClassRes(7:12,1)/250,'b-o');
plot(1:6, regInitClassRes(13:18,1)/250,'k--.');
plot(1:6, regInitClassRes(19:24,1)/250,'m-x');
plot(1:6, regInitClassRes(25:30,1)/250,'-v', 'Color', "#EDB120");
plot(1:6, regInitClassRes(31:36,1)/250,'--', 'Color','#808080');
plot(1:6, regInitClassRes(37:42,1)/250,'->', 'Color', "#D95319");
plot(1:6, regInitClassRes(43:48,1)/250,'-s', 'Color', "#7E2F8E");
plot(1:6, regInitClassRes(49:54,1)/250,'--+', 'Color', "#77AC30");
set(gca, 'xtick', 1:6)
set(gca, 'xticklabel', classNames);
% set(gca, "YTick", 0.9:0.01:1);
% ylim([0.925, 1.005])
ylabel("Robust %");
legend('dropout_G','dropout_H', 'dropout_N', 'jacobian_G','jacobian_H',...
    'jacobian_N', 'L2_G','L2_H', 'L2_N', 'Location','best');
exportgraphics(gca, "plots/comboRes_vs_class.pdf",'ContentType','vector');

% Create figure
figure;
grid; hold on;
plot(1:6, regInitClassRes(1:6,4),'--d', 'Color', "#A2142F");
plot(1:6, regInitClassRes(7:12,4),'b-o');
plot(1:6, regInitClassRes(13:18,4),'k--.');
plot(1:6, regInitClassRes(19:24,4),'m-x');
plot(1:6, regInitClassRes(25:30,4),'-v', 'Color', "#EDB120");
plot(1:6, regInitClassRes(31:36,4),'--', 'Color','#808080');
plot(1:6, regInitClassRes(37:42,4),'->', 'Color', "#D95319");
plot(1:6, regInitClassRes(43:48,4),'-s', 'Color', "#7E2F8E");
plot(1:6, regInitClassRes(49:54,4),'--+', 'Color', "#77AC30");
set(gca, 'xtick', 1:6)
set(gca, 'xticklabel', classNames);
ylabel("Time (s)")
legend('dropout_G','dropout_H', 'dropout_N', 'jacobian_G','jacobian_H',...
    'jacobian_N', 'L2_G','L2_H', 'L2_N', 'Location','best');
exportgraphics(gca, "plots/comboTime_vs_class.pdf",'ContentType','vector');

%% Compare seeds vs classes

seed0 = sum(stats_seed0,1);
seed1 = sum(stats_seed1,1);
seed2 = sum(stats_seed2,1);
seed3 = sum(stats_seed3,1);
seed4 = sum(stats_seed4,1);

% Create figure
figure;
grid;hold on;
plot(1:6, seed0(1,:,1)/450,'r-d');
plot(1:6, seed1(1,:,1)/450,'b-o');
plot(1:6, seed2(1,:,1)/450,'k--');
plot(1:6, seed3(1,:,1)/450,'-x');
plot(1:6, seed4(1,:,1)/450,'-s');
set(gca, 'xtick', 1:6)
set(gca, 'xticklabel', classNames);
% set(gca, "YTick", 0.9:0.01:1);
% ylim([0.925, 1.005])
ylabel("Robust %");
legend('0', '1', '2', '3', '4', 'Location','best');
exportgraphics(gca, "plots/seed_vs_class.pdf",'ContentType','vector');


% Create figure
figure;
grid;hold on;
plot(1:6, seed0(1,:,4)/9,'r-d');
plot(1:6, seed1(1,:,4)/9,'b-o');
plot(1:6, seed2(1,:,4)/9,'k--');
plot(1:6, seed3(1,:,4)/9,'-x');
plot(1:6, seed4(1,:,4)/9,'-s');
set(gca, 'xtick', 1:6)
set(gca, 'xticklabel', classNames);
% set(gca, "YTick", 0.9:0.01:1);
% ylim([0.925, 1.005])
ylabel("Time (s)");
legend('0', '1', '2', '3', '4', 'Location','best');
exportgraphics(gca, "plots/seedTime_vs_class.pdf",'ContentType','vector');

%% Are these models not robuts, or simply harder to verify?
folder_path = 'results_falsify';
results_folder = dir(folder_path);
falseRes = cell(N,2);
for i=3:N+2
    ind_res = load([folder_path, filesep, results_folder(i).name]);
    falseRes{i-2,1} = sum(ind_res.res(:,1)==0);
    falseRes{i-2,2} = results_folder(i).name;
end

%% Helper functions

% Process results per model
function [rob, unk, norob, avg_time] = process_model_res(res, idxs)
    rob   = sum(res(idxs,1) == 1);
    unk   = sum(res(idxs,1) == 2);
    norob = sum(res(idxs,1) == 0);
    avg_time = sum(res(idxs,2))/length(idxs);
end

% Compute average across classes for all models (indexes) in a list
function stats = process_model_classes(res, classes, models_idxs)
    c = size(classes,1);
    m = size(models_idxs,1);
    stats = zeros(m,c,4);
    for i=1:m 
        for j=1:c
            [stats(i,j,1), stats(i,j,2), stats(i,j,3), stats(i,j,4)] = process_model_res(res{models_idxs(i),1}, classes(j,:));
        end
    end
end
