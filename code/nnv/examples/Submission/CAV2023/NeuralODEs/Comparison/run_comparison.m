function run_comparison()

    % Run all nnv reachability experiments
    
    % Set folder path for results
    if ~isfolder('nnvresults')
           mkdir('nnvresults')
    end

    % Turn off figure display
    set(0,'DefaultFigureVisible','off')
    
    %% Run benchmarks 

    cd nnv_reach;
    
    % Cartpole
    CTRNN_Cartpole_reach_short;
    CTRNN_Cartpole_reach_mid;
    CTRNN_Cartpole_reach;

    % Spiral
    spiral2D_linear_reach;
    spiral2D_nonlinear_reach;

    % FPA
    CTRNN_FPA_short;
    CTRNN_FPA_mid;
    CTRNN_FPA_reach;

    %% Visualize results

    cd ..;
    plot_spiralL(3);
    plot_spiralNL(3);
    process_results_neuralode;

end
