function plot_spiralL(scenario)

    %% Plot one scenario for the linear spiral 2D
    % scenario = 1, 2 or 3, equivalent to 0.01, 0.05 or 0.1    
    
    %% NNV
    if scenario == 3
        sp_nnv = load("nnvresults/spiral_0.1.mat");
    elseif scenario == 2
        sp_nnv = load("nnvresults/spiral_0.05.mat");
    elseif scenario == 1
        sp_nnv = load("nnvresults/spiral_0.01.mat");
    else
        error("Wrong scenario specified. Scenario must be 1, 2, or 3")
    end
    
    
    %% Generate Plots
    
    figure;
    hold on;
    Star.plotBoxes_2D_noFill(sp_nnv.Rall(1:end-1),1,2,'b'); % NNV 
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');
    legend('off')
    grid;
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 14; % Set font size of axis
    ax.YAxis.FontSize = 14;
    legend('off')
    hold on;
    pj = plot(2, 0,'b');
    legend(pj,{'NNV 2.0'},"Location","best",'FontSize',13);
    
    % Save figure
    if is_codeocean
        saveas(gcf,"/results/logs/spiralL_"+string(scenario)+".pdf");
    else
        saveas(gcf,"spiralL_"+string(scenario)+".pdf");
    end

end