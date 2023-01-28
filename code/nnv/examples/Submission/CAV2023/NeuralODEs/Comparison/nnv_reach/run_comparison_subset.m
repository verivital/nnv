function run_comparison_subset()

    %% Set folder path for results
    if ~exist('../nnvresults', 'dir')
           mkdir('../nnvresults')
    end
    
    %% Run benchmarks tool compare
    
    % Spiral
    spiral2D_linear_reach;
    
    %% FPA
    fpa_reach_short;
    fpa_reach_mid;
    fpa_reach_long;
    
end
