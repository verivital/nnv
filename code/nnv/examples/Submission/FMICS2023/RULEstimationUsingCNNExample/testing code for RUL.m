load('TEDS_results.mat', 'X')
load('TEDS_results.mat', 'Y')
load('TEDS_results.mat', 'Y_pred_52')
load('TEDS_results.mat', 'dataset_TEDS')
load('TEDS_results.mat', 'net_TEDS')
load('TEDS_results.mat', 'perturbations')
load('TEDS_results.mat', 'feature_to_attack')

nnvnet_TEDS = matlab2nnv(net_TEDS);

%% i represenrts the dataset for different sequence lengths
% i == 1 -> seq length 10
% i == 2 -> seq length 20
% i == 3 -> seq length 30

i = 1;

for j = 1: length(perturbations)
    
    [LB_n1{i,j},UB_n1{i,j},T_n1{i,j},PR_n1{i,j},POR_n1{i,j},T_avg_n1{i,j},T_sum_n1{i,j}] = reachabilityNoise(nnvnet,dataset_cons,i,perturbations(j),'SFSI',feature_to_attack);
   % [LB_n2{i,j},UB_n2{i,j},T_n2{i,j},PR_n2{i,j},POR_n2{i,j},T_avg_n2{i,j},T_sum_n2{i,j}] = reachabilityNoise(nnvnet,dataset_cons,i,perturbations(j),'SFAI',feature_to_attack);
   % [LB_n3{i,j},UB_n3{i,j},T_n3{i,j},PR_n3{i,j},POR_n3{i,j},T_avg_n3{i,j},T_sum_n3{i,j}] = reachabilityNoise(nnvnet,dataset_cons,i,perturbations(j),'MFSI',feature_to_attack);

end
 