%% Run all benchmarks with NNV
% Random examples
cd RandomEx
run run_all.m
cd ..;

% MNIST
cd MNIST
run run_mnist.m
cd ..

% ACC
cd ACC
run verify_acc_orig.m
run verify_acc_plant.m
run verify_acc_tanhplant.m
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