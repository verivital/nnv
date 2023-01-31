function run_neuralode_subset()

    %% Run all benchmarks with NNV

    % Random examples
    cd RandomEx
    reach_XS;
    reach_S;
    reach_M;
    randomEx_table_subset;
    cd ..;
    
    % MNIST
    cd MNIST
    cd ffnn
    run_ffnn_all;
    run_ffnn_all_inf;
    cd ..
    mnist_table_subset;
    cd ..
    
    % ACC
    cd ACC
    reach_acc_orig;
    reach_acc_plant;
    reach_acc_tanhplant;
    cd ..
    
    % NODE Comparison
    cd Comparison
    cd nnv_reach
    run_comparison_subset
    cd ..
    tool_comparison_table_subset.m % create Table 3
    process_results_neuralode_subset.m % plot FPA
    plot_spiralL;
    cd ..;

end