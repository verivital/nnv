function run_neuralode()

    % Run all neuralode benchmarks with NNV
    
    % Random examples
    cd RandomEx
    run_neuralode_random;
    cd ..;
    
    % MNIST
    cd MNIST
    run run_mnist.m
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
    run run_nnv.m
    cd ..
    run generate_table.m % create Table 3
    run process_results.m % plot FPA, Cartpole and Damped Oscillator
    run plot_spiralL.m % plot spiral results
    run plot_spiralNL.m % plot nonlinear spiral results
    cd ../..

end