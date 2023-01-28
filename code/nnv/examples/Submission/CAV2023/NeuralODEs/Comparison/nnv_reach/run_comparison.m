function run_comparison()

    % Run all nnv reachability experiments
    %% Set folder path for results
    if ~exist('../nnvresults', 'dir')
           mkdir('../nnvresults')
    end
%     addpath('benchmark_dynamics/');
    % Turn off figure display
    set(0,'DefaultFigureVisible','off')
    
    %% Run benchmarks tool compare
    
    % Cartpole
    run CTRNN_Cartpole_reach_short.m
    run CTRNN_Cartpole_reach_mid.m
    run CTRNN_Cartpole_reach.m
    % Spiral
    run spiral2D_linear_reach.m
    run spiral2D_nonlinear_reach.m
    % FPA
    run CTRNN_FPA_short.m
    run CTRNN_FPA_mid.m
    run CTRNN_FPA_reach.m
    
    % Damped Oscillator
    cd DampedOscillator/
    run_compare_dim()
    cd ..

end
