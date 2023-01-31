function run_neuralode()

    % Run all neuralode benchmarks with NNV
    
    % Random examples
    cd RandomEx
    run_neuralode_random;
    cd ..;
    
    % MNIST
    cd MNIST
    run_mnist;
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
    run_comparison;

    cd ..
    tool_comparison_table.m    % create Table 3
    process_results_neuralode; % plot FPA, Cartpole and Damped Oscillator
    plot_spiralL;              % plot spiral results
    plot_spiralNL;             % plot nonlinear spiral results
    cd ..;

end