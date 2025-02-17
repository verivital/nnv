%% Verify all possible 2D classification models for medmnist data

medmnist_path = "data/organcmnist.mat"; % path to data
model_path = "models/model_organcmnist.mat"; % path to trained models

disp("Begin verification of organCMNIST");

% Load data
load(medmnist_path);

% data to verify (test set)
test_images = permute(test_images, [2 3 4 1]);
test_labels = test_labels + 1;

% load network
load(model_path);
net = matlab2nnv(net);

% adversarial attack
adv_attack = struct;
adv_attack.Name = "linf";
epsilons = [1,2,3,4];

for ep = epsilons
    adv_attack.epsilon = ep; % {epsilon} color values
    
    % select images to verify
    % N = 50;
    N = 10;
    inputs = test_images(:,:,:,1:N);
    targets = test_labels(1:N);
    
    % verify images
    results = verifyDataset(net, inputs, targets, adv_attack);
    
    % save results
    save("results/verification_organcmnist_mini_"+string(ep)+".mat", "results", "adv_attack");

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
    
    disp("|========================================================================|")
    disp(' ');

end

%% Create table results

results_1 = load("results/verification_organcmnist_mini_1.mat");
results_2 = load("results/verification_organcmnist_mini_2.mat");
results_3 = load("results/verification_organcmnist_mini_3.mat");
results_4 = load("results/verification_organcmnist_mini_4.mat");

% Robust
r1_rb = sum(results_1.results(1,:) == 1)/size(results_1.results,2);
r2_rb = sum(results_2.results(1,:) == 1)/size(results_2.results,2);
r3_rb = sum(results_3.results(1,:) == 1)/size(results_3.results,2);
r4_rb = sum(results_4.results(1,:) == 1)/size(results_4.results,2);

% Not Robust
r1_nrb = sum(results_1.results(1,:) == 0)/size(results_1.results,2);
r2_nrb = sum(results_2.results(1,:) == 0)/size(results_2.results,2);
r3_nrb = sum(results_3.results(1,:) == 0)/size(results_3.results,2);
r4_nrb = sum(results_4.results(1,:) == 0)/size(results_4.results,2);

% Unknown
r1_unk = sum(results_1.results(1,:) == 2)/size(results_1.results,2);
r2_unk = sum(results_2.results(1,:) == 2)/size(results_2.results,2);
r3_unk = sum(results_3.results(1,:) == 2)/size(results_3.results,2);
r4_unk = sum(results_4.results(1,:) == 2)/size(results_4.results,2);

% Misclassified
r1_ms = sum(results_1.results(1,:) == -1)/size(results_1.results,2);
r2_ms = sum(results_2.results(1,:) == -1)/size(results_2.results,2);
r3_ms = sum(results_3.results(1,:) == -1)/size(results_3.results,2);
r4_ms = sum(results_4.results(1,:) == -1)/size(results_4.results,2);

% Average Time
r1_time = sum(results_1.results(2,:))/size(results_1.results,2);
r2_time = sum(results_2.results(2,:))/size(results_2.results,2);
r3_time = sum(results_3.results(2,:))/size(results_3.results,2);
r4_time = sum(results_4.results(2,:))/size(results_4.results,2);

resTable = {[1;2;3;4], [r1_rb;r2_rb;r3_rb;r4_rb],...
    [r1_nrb;r2_nrb;r3_nrb;r4_nrb],...
    [r1_unk;r2_unk;r3_unk;r4_unk],...
    [r1_ms;r2_ms;r3_ms;r4_ms],...
    [r1_time;r2_time;r3_time;r4_time]};

VarNames = ["Epsilon", "Robust", "Not Robust", "Unknown", "Misclass.", "Avg Time (s)"];

%%
f = figure;
uit = uitable(f, 'Data', [resTable{1}, resTable{2}, resTable{3}, resTable{4}, resTable{5}, resTable{6}]);
uit.ColumnName = VarNames;
uit.RowName=[]; %removing default row numbering as in your uitable
table_extent = get(uit,'Extent');
set(uit,'Position',[1 1 table_extent(3) table_extent(4)])
figure_size = get(f,'outerposition');
desired_fig_size = [figure_size(1)-15 figure_size(2) table_extent(3)+55 table_extent(4)+95];
set(f,'outerposition', desired_fig_size);
print('-dpng', 'results_table.png', '-r0'); %-r0 prints at screen resolution

