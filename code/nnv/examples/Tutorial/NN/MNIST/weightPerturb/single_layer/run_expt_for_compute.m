clear
empty_results_file = 'results/MNIST_MLP_empty_results_file.yaml';
results_file = 'results/MNIST_MLP';
copyfile(empty_results_file, [results_file '.yaml'])
do_not_clear = 1;
n_layers_to_run_for_from_yaml_file = 3;
expt_path = results_file;
conv_expt_any_layer