function plot_spiralL(scenario)

    %% Plot one scenario for the linear spiral 2D
    % scenario = 1, 2 or 3, equivalent to 0.01, 0.05 or 0.1
    
    %% Gotube (no results)
    
    
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
    
    %% JuliaReach
    if is_codeocean
        julia_path = "/results/logs/juliaresults/";
    else
        julia_path = "juliareach/results/";
    end
    sp_julia = load(julia_path + "spiral" + string(scenario)+ ".mat");
    
    
    %% Generate Plots
    
    figure;
    hold on;
    Star.plotBoxes_2D_noFill(sp_nnv.Rall(1:end-1),1,2,'k'); % NNV 
    plot_juliareach(sp_julia,1,2,'b'); % Juliareach
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');
    legend('off')
    grid;
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 14; % Set font size of axis
    ax.YAxis.FontSize = 14;
    % Flowstar reach sets
    plot_flowstar(scenario)
    legend('off')
    hold on;
    pn = plot(2, 0,'k');
    pj = plot(2, 0,'b');
    pf = plot(2, 0,'Color',[0 0.4 0]);
    legend([pf pn pj],{'Flow*', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
    
    % Save figure
    if is_codeocean
        saveas(gcf,"/results/logs/spiralL_compare_"+string(scenario)+".pdf");
    else
        saveas(gcf,"spiralL_compare_"+string(scenario)+".pdf");
    end
    pause(0.01); % Add a pause to ensure figure is saved
    xlim([0.23 0.49]);
    ylim([0.48 0.7]);
    if is_codeocean
        saveas(gcf,"/results/logs/spiralL_compare_"+string(scenario)+"_zoom.pdf");
    else
        saveas(gcf,"spiralL_compare_"+string(scenario)+"_zoom.pdf");
    end

end

%% Helper Functions

function plot_flowstar(scenario)
    if is_codeocean
        if scenario == 3
            run /results/logs/flowresults/outputs/spiralL3.m
        elseif scenario == 2
            run /results/logs/flowresults/outputs/spiralL2.m
        else
            run /results/logs/flowresults/outputs/spiralL1.m
        end
    else
        if scenario == 3
            run flowresults/results/spiralL3.m
        elseif scenario == 2
            run flowresults/results/spiralL2.m
        else
            run flowresults/results/spiralL1.m
        end
    end
end
