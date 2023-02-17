%% Run all benchmarks with NNV
% Random examples
cd RandomEx
run reach_XS.m
run reach_S.m
run reach_M.m
run generate_table_subset.m
cd ..;

% MNIST
cd MNIST
cd ffnn
run run_all.m
run run_all_inf
cd ..
run create_table_subset.m
cd ..

% ACC
cd ACC
run verify_acc_orig.m
run verify_acc_plant.m
cd ..

% NODE Comparison
cd Comparison
cd nnv_reach
run run_subset.m
cd ..
run generate_table_subset.m % create Table 3
run process_results_subset.m % plot FPA
run plot_spiralL.m
cd ../..