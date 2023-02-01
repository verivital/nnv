function run_neuralode()

    % Run all neuralode benchmarks with NNV
    
    % Random examples
    cd RandomEx
    if ~exist('results', 'dir')
       mkdir('results')
    end
    reach_M;
    reach_XXL;
    cd ..;
    
    % MNIST
    cd MNIST
    run_mnist_inf;
    mnist_table;
    cd ..
    
    % ACC
    cd ACC
    reach_acc_orig;
    reach_acc_plant;
    reach_acc_tanhplant;
    cd ..
    
    % NODE Comparison
    cd Comparison
    run_comparison;
    cd ..;

end