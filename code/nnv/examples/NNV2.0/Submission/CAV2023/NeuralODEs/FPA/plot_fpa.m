function plot_fpa()

    %% Process and plot all FPA and Cartpole reachability results
    
    %% 3 - Load NNV results

    % Damped Forced Pendulum
    fpa_nnv = "fpa_reach.mat";
    fpa_nnv = load(fpa_nnv);
    
    
    %% 5 - Plot FPA reach sets
    f_fpa = figure;
    hold on;
    Star.plotBoxes_2D_noFill(fpa_nnv.Rall,1,2,'b');
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');
    % Add initial state to generate legends x0 = [0.21535, -0.58587]
    legend('off')
    grid;
    pj = plot(0.21535, -0.58587,'b');
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 15; % Set font size of axis
    ax.YAxis.FontSize = 15;
    legend(pj,{'NNV'},"Location","best",'FontSize',14);
    if is_codeocean
        exportgraphics(f_fpa,'/results/logs/fpa.pdf', 'ContentType', 'vector');
        saveas(f_fpa,'/results/logs/fpa.png');
    else
        exportgraphics(f_fpa,'../../fpa.pdf','ContentType', 'vector');
    end

end