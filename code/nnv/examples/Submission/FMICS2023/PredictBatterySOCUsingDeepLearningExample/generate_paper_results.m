%% this file generates the outputs presented in the paper:

dataset_SOC = createSOCDataset(X, Y, Y_Pred_n10degC, [10, 20, 30], 5, "cons");
nnvnet_SOC = matlab2nnv(net);
noise_extent_SOC = [0.01, 0.025, 0.05,0.1];
feature_to_attack = 3; % can be random

for i = 1 : length(dataset_SOC.input)
    for j = 1 : length(noise_extent_SOC)
        [LB_MFAI{i,j},UB_MFAI{i,j},T_MFAI{i,j},PR_MFAI{i,j},POR_MFAI{i,j},T_avg_MFAI{i,j},T_sum_MFAI{i,j}] = reachabilityNoiseForSOC(nnvnet_SOC,dataset_SOC,i,noise_extent_SOC(j),'MFAI');
        [LB_MFSI{i,j},UB_MFSI{i,j},T_MFSI{i,j},PR_MFSI{i,j},POR_MFSI{i,j},T_avg_MFSI{i,j},T_sum_MFSI{i,j}] = reachabilityNoiseForSOC(nnvnet_SOC,dataset_SOC,i,noise_extent_SOC(j),'MFSI');
        [LB_SFAI{i,j},UB_SFAI{i,j},T_SFAI{i,j},PR_SFAI{i,j},POR_SFAI{i,j},T_avg_SFAI{i,j},T_sum_SFAI{i,j}] = reachabilityNoiseForSOC(nnvnet_SOC,dataset_SOC,i,noise_extent_SOC(j),'SFAI',feature_to_attack);
        [LB_SFSI{i,j},UB_SFSI{i,j},T_SFSI{i,j},PR_SFSI{i,j},POR_SFSI{i,j},T_avg_SFSI{i,j},T_sum_SFSI{i,j}] = reachabilityNoiseForSOC(nnvnet_SOC,dataset_SOC,i,noise_extent_SOC(j),'SFSI',feature_to_attack);
    end
end