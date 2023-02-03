function plot_spiralNL(scenario)
    
    %% Plot one scenario for the nonlinear spiral 2D
    
    %% NNV
    if scenario == 1
        sp_nnv = load("nnvresults/spiral_nl_0.01.mat");
    elseif scenario == 2
        sp_nnv = load("nnvresults/spiral_nl_0.05.mat");
    else
        sp_nnv = load("nnvresults/spiral_nl_0.1.mat");
    end
    
    %% Generate Plots
    
    % Create figure
    f1 = figure;
    hold on;
    Star.plotBoxes_2D_noFill(sp_nnv.Rall(1:end-1),1,2,'b');
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');
    legend('off')
    grid;
    pj = plot(2, 0,'b');
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 14; % Set font size of axis
    ax.YAxis.FontSize = 14;
    legend(pj,{'NNV'},"Location","best",'FontSize',13);
    if is_codeocean
        saveas(f1,"/results/logs/spiral_"+string(scenario)+"_nl.pdf");
    else
        saveas(f1,"spiral_"+string(scenario)+"_nl.pdf");
    end

end