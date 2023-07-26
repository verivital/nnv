%% this file genertes the Fig. 2 captioned 
% "Allowable (blue) and reachable (red) bounds for battery SOC dataset for
% 100 consecutive time steps and 2 different SFAI noise values 1% (upper), and
% 2.5% (lower) respectively" 
load SOC_results.mat;

for i = 1:3
    for j= 1:100
        LB{i,1}(j)=dataset_SOC.allowableLB{i, 1}{j, 1}(end);
        UB{i,1}(j)=dataset_SOC.allowableUB{i, 1}{j, 1}(end);
    end
end
y_actual = (LB{1,1}+UB{1,1})/2;
err1 = y_actual- LB{1,1};

figure
subplot(3,1,1)
errorbar(x,y_actual,err1,"MarkerSize",10,"LineStyle","none", "LineWidth",2)
hold on
err2 = (UB_SFAI{3,1} - LB_SFAI{3,1})/2;
y_pred = (UB_SFAI{3,1} + LB_SFAI{3,1})/2;
errorbar(x,y_pred,err2,err2,"MarkerSize",10,"LineStyle","none","LineWidth",2)
legend('Allowable Bounds','Reachable Bounds')   

subplot(3,1,2)
errorbar(x,y_actual,err1,"MarkerSize",10,"LineStyle","none","LineWidth",2)
hold on
err2 = (UB_SFAI{3,2} - LB_SFAI{3,2})/2;
y_pred = (UB_SFAI{3,2} + LB_SFAI{3,2})/2;
errorbar(x,y_pred,err2,err2,"MarkerSize",10,"LineStyle","none","LineWidth",2)
legend('Allowable Bounds','Reachable Bounds')

subplot(3,1,3)
errorbar(x,y_actual,err1,"MarkerSize",10,"LineStyle","none","LineWidth",2)
hold on
err2 = (UB_SFAI{3,3} - LB_SFAI{3,3})/2;
y_pred = (UB_SFAI{3,3} + LB_SFAI{3,3})/2;
errorbar(x,y_pred,err2,err2,"MarkerSize",10,"LineStyle","none","LineWidth",2)
legend('Allowable Bounds','Reachable Bounds')