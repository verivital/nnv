function run_mnist()
    
    %% Run all MNIST related experiments
    cd cnn;
    run_cnn_all;
    run_cnn_all_inf;
    eval_cnn_all;
        
    cd ..;
    cd ffnn;
    run_ffnn_all;
    run_ffnn_all_inf;
    eval_ffnn_all;
        
    cd ..
    mnist_table;

end