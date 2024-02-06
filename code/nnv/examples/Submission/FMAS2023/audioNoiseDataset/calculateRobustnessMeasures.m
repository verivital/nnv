nnvnet = matlab2nnv(net);
load('data.mat',wNoiseTest_ex);
load('data.mat',bNoiseTest_ex);
load('data.mat',pNoiseTest_ex);

for i = 1: length(percent)
for j = 1:1
[robust_SFSI{i,j},T_SFSI{i,j},PR_SFSI{i,j},T_avg_SFSI{i,j},T_sum_SFSI{i,j}] = adversarialReachability(nnvnet,wNoiseTest_ex,1,percent(i),"SFSI",j);
[robust_SFAI{i,j},T_SFAI{i,j},PR_SFAI{i,j},T_avg_SFAI{i,j},T_sum_SFAI{i,j}] = adversarialReachability(nnvnet,wNoiseTest_ex,1,percent(i),"SFAI",j);
[robust_MFSI{i,j},T_MFSI{i,j},PR_MFSI{i,j},T_avg_MFSI{i,j},T_sum_MFSI{i,j}] = adversarialReachability(nnvnet,wNoiseTest_ex,1,percent(i),"MFSI",j);
[robust_MFAI{i,j},T_MFAI{i,j},PR_MFAI{i,j},T_avg_MFAI{i,j},T_sum_MFAI{i,j}] = adversarialReachability(nnvnet,wNoiseTest_ex,1,percent(i),"MFAI",j);
end
end
for i = 1: length(percent)
for j = 2:2
[robust_SFSI{i,j},T_SFSI{i,j},PR_SFSI{i,j},T_avg_SFSI{i,j},T_sum_SFSI{i,j}] = adversarialReachability(nnvnet,bNoiseTest_ex,1,percent(i),"SFSI",j);
[robust_SFAI{i,j},T_SFAI{i,j},PR_SFAI{i,j},T_avg_SFAI{i,j},T_sum_SFAI{i,j}] = adversarialReachability(nnvnet,bNoiseTest_ex,1,percent(i),"SFAI",j);
[robust_MFSI{i,j},T_MFSI{i,j},PR_MFSI{i,j},T_avg_MFSI{i,j},T_sum_MFSI{i,j}] = adversarialReachability(nnvnet,bNoiseTest_ex,1,percent(i),"MFSI",j);
[robust_MFAI{i,j},T_MFAI{i,j},PR_MFAI{i,j},T_avg_MFAI{i,j},T_sum_MFAI{i,j}] = adversarialReachability(nnvnet,bNoiseTest_ex,1,percent(i),"MFAI",j);
end
end
for i = 1: length(percent)
for j = 3:3
[robust_SFSI{i,j},T_SFSI{i,j},PR_SFSI{i,j},T_avg_SFSI{i,j},T_sum_SFSI{i,j}] = adversarialReachability(nnvnet,pNoiseTest_ex,1,percent(i),"SFSI",j);
[robust_SFAI{i,j},T_SFAI{i,j},PR_SFAI{i,j},T_avg_SFAI{i,j},T_sum_SFAI{i,j}] = adversarialReachability(nnvnet,pNoiseTest_ex,1,percent(i),"SFAI",j);
[robust_MFSI{i,j},T_MFSI{i,j},PR_MFSI{i,j},T_avg_MFSI{i,j},T_sum_MFSI{i,j}] = adversarialReachability(nnvnet,pNoiseTest_ex,1,percent(i),"MFSI",j);
[robust_MFAI{i,j},T_MFAI{i,j},PR_MFAI{i,j},T_avg_MFAI{i,j},T_sum_MFAI{i,j}] = adversarialReachability(nnvnet,pNoiseTest_ex,1,percent(i),"MFAI",j);
end
end

GPR_MFAI = sum(cell2mat(PR_MFAI)')/3;
GPR_MFSI = sum(cell2mat(PR_MFSI)')/3;
GPR_SFSI = sum(cell2mat(PR_SFSI)')/3;
GPR_SFAI = sum(cell2mat(PR_SFAI)')/3;
GTavg_MFAI = sum(cell2mat(T_avg_MFAI)')/3;
GTavg_MFSI = sum(cell2mat(T_avg_MFSI)')/3;
GTavg_SFSI = sum(cell2mat(T_avg_SFSI)')/3;
GTavg_SFAI = sum(cell2mat(T_avg_SFAI)')/3;
GTsum_MFAI = sum(cell2mat(T_sum_MFAI)')/3;
GTsum_MFSI = sum(cell2mat(T_sum_MFSI)')/3;
GTsum_SFSI = sum(cell2mat(T_sum_SFSI)')/3;
GTsum_SFAI = sum(cell2mat(T_sum_SFAI)')/3;