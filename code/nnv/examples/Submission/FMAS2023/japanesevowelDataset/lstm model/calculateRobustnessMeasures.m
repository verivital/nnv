percent = [0.5, 0.6, 0.7, 0.8, 0.9];
nnvnet = matlab2nnv(net);
YTest = double(YTest_correct);
for j = 1: length(percent)
    for i = 1:length(XTest_correct)
        input = XTest_correct{i,1};
        classIndex = YTest(i);
        [robust_SFSI(j,i),T_SFSI(j,i)] = adversarialReachability(nnvnet,input,1,percent(j),"SFSI",classIndex);
        [robust_SFAI(j,i),T_SFAI(j,i)] = adversarialReachability(nnvnet,input,1,percent(j),"SFAI",classIndex);
        [robust_MFSI(j,i),T_MFSI(j,i)] = adversarialReachability(nnvnet,input,1,percent(j),"MFSI",classIndex);
        [robust_MFAI(j,i),T_MFAI(j,i)] = adversarialReachability(nnvnet,input,1,percent(j),"MFAI",classIndex);
    end
    PR_SFSI(1,j) = sum(robust_SFSI(j,:)==1)/length(XTest_correct);
    T_sum_SFSI(1,j) = sum(T_SFSI(j,:));
    PR_SFAI(1,j) = sum(robust_SFAI(j,:)==1)/length(XTest_correct);
    T_sum_SFAI(1,j) = sum(T_SFAI(j,:));
    PR_MFSI(1,j) = sum(robust_MFSI(j,:)==1)/length(XTest_correct);
    T_sum_MFSI(1,j) = sum(T_MFSI(j,:));
    PR_MFAI(1,j) = sum(robust_MFAI(j,:)==1)/length(XTest_correct);
    T_sum_MFAI(1,j) = sum(T_MFAI(j,:));
end