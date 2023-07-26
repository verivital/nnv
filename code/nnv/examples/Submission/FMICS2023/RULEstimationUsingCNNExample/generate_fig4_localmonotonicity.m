%% this file generates fig.4 (local monotonicity), of the paper 
load ('TEDS_results.mat','X','Y','Y_pred_52');

plot_localmonotonicity(Y,Y_pred_52);

function plot_localmonotonicity(true_RUL, predicted_RUL)
%% this function 'plot_localmonotonicity' plots the true RUL, the predicted
% RUL and the linearized version of the predicted RUL for n time-instances 

    predicted_lin = fitLinearLine(true_RUL);

    figure;
    
    x = 1 : length(true_RUL);
    plot(x, true_RUL,'b');
    hold on
    plot(x, predicted_RUL, 'r');
    hold on
    plot(x, predicted_lin, '-');
    legend('True RUL','Predicted RUL','Linearized Predicted RUL')
    xlabel('Time stamp (Test data sequence)');
    ylabel('RUL (Cycles)');
    title('RUL for Test engine #52 (52 case)');
    set(gca, 'FontSize', 22);
    
end

function y_lin = fitLinearLine(y)
    x = 1:length(y);
    c = polyfit(x,y,1);
    y_lin = polyval(c,x);
end