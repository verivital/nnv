function process_results_neuralode()

    %% Process and plot all FPA and Cartpole reachability results
    
    %% 3 - Load NNV results
    nnv_path = "nnvresults/";
    % Damped Forced Pendulum
    fpa_nnv = "fpa_reach.mat";
    fpa_nnv = load(nnv_path+fpa_nnv);
    % Cartpole
    cartpole_nnv = "cartpole_reach.mat";
    cartpole_nnv = load(nnv_path+cartpole_nnv);
    
    
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
    legend(pj,{'NNV 2.0'},"Location","best",'FontSize',14);
    if is_codeocean
        saveas(f_fpa,'/results/logs/fpa.pdf');
    else
        saveas(f_fpa,'fpa.pdf');
    end
    
    %% 6 - Plot Cartpole reachsets
    f_cp = figure;
    hold on;
    Star.plotBoxes_2D_noFill(cartpole_nnv.Rall,1,2,'b');
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');
    legend('off');
    grid;
    pj = plot(0,0,'b');
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 15; % Set font size of axis
    ax.YAxis.FontSize = 15;
    legend(pj,{'NNV 2.0'},"Location","best",'FontSize',14);
    if is_codeocean
        saveas(f_fpa,'/results/logs/cartpole.pdf');
    else
        saveas(f_cp,'cartpole.pdf');
    end
    

end