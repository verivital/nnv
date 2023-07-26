function plot_localmonotonicity(LB, UB)
%% this function 'plot_localmonotonicity' plots the true RUL, the predicted
% RUL and the linearized version of the predicted RUL for n time-instances 

    LB_lin = fitLinearLine(LB);
    UB_lim = fitLinearLine(UB);

    figure;
    
end

function y_lin = fitLinearLine(y)
    x = 1:length(y);
    c = polyfit(x,y,1);
    y_lin = polyval(c,x);
end