function run_neuralode()

    % Run all neuralode benchmarks with NNV

    disp("Running MNIST experiments...");
    
    % MNIST
    cd MNIST
    run_mnist_inf;
    cd ..
    
    disp("Running ACC example...");

    % ACC
    cd ACC
    reach_acc_nonlinear;
    cd ..

    disp("Running FPA example...")
    
    % FPA
    cd FPA;
    fpa_reach;
    plot_fpa;
    cd ..;

end