function run_neuralode()

    % Run all neuralode benchmarks with NNV
    
    % MNIST
    cd MNIST
    run_mnist_inf;
    mnist_table;
    cd ..
    
    % ACC
    cd ACC
    reach_acc_nonlinear;
    cd ..
    
    % NODE Comparison
    cd FPA;
    fpa_reach;
    plot_fpa;
    cd ..;

end