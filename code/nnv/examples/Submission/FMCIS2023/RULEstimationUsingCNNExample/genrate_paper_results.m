

% dataset_cons = createDataset(X, Y, Y_pred_52, [10, 20, 30], "cons");
% nnvnet = matlab2nnv(net);
perturbations = [0.1, 0.25, 0.5];
% SFSI
for i = 3:length(dataset_cons.input)
    for j = 1: length(perturbations)
    %[LB_n1{i,j},UB_n1{i,j},T_n1{i,j},PR_n1{i,j},POR_n1{i,j},T_avg_n1{i,j},T_sum_n1{i,j}] = reachabilityNoise(nnvnet,dataset_cons,i,perturbations(j),'SFSI',3);
    [LB_n2{i,j},UB_n2{i,j},T_n2{i,j},PR_n2{i,j},POR_n2{i,j},T_avg_n2{i,j},T_sum_n2{i,j}] = reachabilityNoise(nnvnet,dataset_cons,i,perturbations(j),'SFAI',3);
%     [LB_n3{i,j},UB_n3{i,j},T_n3{i,j},PR_n3{i,j},POR_n3{i,j},T_avg_n3{i,j},T_sum_n3{i,j}] = reachabilityNoise(nnvnet,dataset_cons,i,perturbations(j),'MFSI',3);

    end
end