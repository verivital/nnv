function plot_spiralNL(scenario)
    %% Plot one scenario for the nonlinear spiral 2D
    % Flowstar timed out, no results
    
    %% Gotube
    if is_codeocean
        gotube_path = "/results/logs/tuberesults/";
    else
        gotube_path = "GoTube/saved_outputs/";
    end
    
    % get file name
    if scenario == 1
        sp_gt = "spiralNL_10.0_0.01_1000_0.01_0.01_1.5_GoTube.txt";
    elseif scenario == 2
        sp_gt = "spiralNL_10.0_0.01_1000_0.05_0.01_1.5_GoTube.txt";
    elseif sceario == 3
        sp_gt = "spiralNL_10.0_0.01_1000_0.1_0.01_1.5_GoTube.txt";
    else
        error("Wrong scenario specified. Scenario must be 1, 2, or 3")
    end
    % load data
    sp_gt = gotube_path+sp_gt; 
    fileID = fopen(sp_gt,'r');
    sp_gt = fscanf(fileID, "%f", [8 1001]);
    fclose(fileID);
    
    %% NNV
    if scenario == 1
        sp_nnv = load("results/spiral_nl_0.01.mat");
    elseif scenario == 2
        sp_nnv = load("results/spiral_nl_0.05.mat");
    else
        sp_nnv = load("nnvresults/spiral_nl_0.1.mat");
    end
    
    %% JuliaReach
    if is_codeocean
        julia_path = "/results/logs/juliaresults/";
    else
        julia_path = "juliareach/results/";
    end
    
    sp_julia = load(julia_path + "spiral" + string(scenario)+".mat");
    
    
    %% Generate Plots
    
    % Create figure
    f1 = figure;
    hold on;
    plot_gotube(sp_gt,1,2,2,'r'); 
    Star.plotBoxes_2D_noFill(sp_nnv.Rall(1:end-1),1,2,'k');
    plot_juliareach(sp_julia,1,2,'b');
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');
    legend('off')
    grid;
    pg = plot(2, 0,'r');
    pn = plot(2, 0,'k');
    pj = plot(2, 0,'b');
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 14; % Set font size of axis
    ax.YAxis.FontSize = 14;
    legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
    if is_codeocean
        saveas(f1,"/results/logs/spiral_compare_"+string(scenario)+"_nl.pdf");
    else
        saveas(f1,"spiral_compare_"+string(scenario)+"_nl.pdf");
    end
    
    % Create figure
    f1 = figure;
    hold on;
    plot_gotube(sp_gt(:,end),1,2,2,'r'); 
    Star.plotBoxes_2D_noFill(sp3_nnv.Rall(end-1),1,2,'k');
    plot_juliareach(sp_julia,1,2,'b',"last");
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');
    legend('off')
    grid;
    pg = plot(2, 0,'r');
    pn = plot(2, 0,'k');
    pj = plot(2, 0,'b');

    ax = gca; % Get current axis
    ax.XAxis.FontSize = 14; % Set font size of axis
    ax.YAxis.FontSize = 14;
    set_axis_limits(scenario);
    legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
    if is_codeocean
        saveas(f1,"/results/logs/spiral_compare_"+string(scenario)+"_nl_last.pdf");
    else
        saveas(f1,"spiral_compare_"+string(scenario)+"_nl_last.pdf");
    end

end

%% Helper functions

function set_axis_limits(scenario)
    switch scenario
        case 1
            xlim([0.52 0.62])
            ylim([0.46 0.51])
        case 2
            xlim([0.45 0.7])
            ylim([0.41 0.56])
        case 3
            xlim([0.2 0.9])
            ylim([0.2 0.8])
        otherwise
            error("Wrong scenario specified. Scenario must be 1, 2, or 3")
    end
end