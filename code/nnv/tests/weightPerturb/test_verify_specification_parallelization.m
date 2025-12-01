
% the result should be the same using the parallelized as well as the
% non-parallelized version of verify_specification

clear

nn_name = 'Small_ConvNet';
folder = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Submission', filesep, 'CAV2020_ImageStar', filesep, 'MNIST_NETS', filesep, 'Architecture', filesep];
mnist_model = load([folder nn_name '.mat']);
matlabnet = mnist_model.net;

% Create NNV model
net = matlab2nnv(matlabnet);
netbs = getByteStreamFromArray(net);    % bytestream to copy this network
    % for weight perturbations in parfor loop iterations
        
% Define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approximate method
reachOptions.numCores = 1;
reachOptions.device = 'cpu';
reachOptions.free_mem_frac_for_LPs_in_verify_specification = 0.1;

reachOptions_parallel = reachOptions;
reachOptions_parallel.layer_specific_numCores = dictionary( ...
    ["start", "Conv2DLayer", "ReluLayer", "MaxPooling2DLayer", "end"], ...
    [  10   ,      10      ,     10     ,         10         ,  10  ]);
reachOptions_parallel.dis_opt = 'display';

[images, labels] = load_images(database = "mnist", ...
    n = 10, ...
    matlabnet = matlabnet);
img = images{1};

layer_ind = 13;

%% test 1: small perturbation, robust case
frac_of_weight_range = 0.001;
test_for_given_frac_of_weight_range(netbs, frac_of_weight_range, layer_ind, img, reachOptions, reachOptions_parallel);

%% test 2: large perturbation, non-robust case
frac_of_weight_range = 0.01;
test_for_given_frac_of_weight_range(netbs, frac_of_weight_range, layer_ind, img, reachOptions, reachOptions_parallel);

delete(gcp('nocreate'));

%% test code
function test_for_given_frac_of_weight_range(netbs, frac_of_weight_range, layer_ind, img, reachOptions, reachOptions_parallel)
    delete(gcp('nocreate'));
    net = getArrayFromByteStream(netbs);
    layer_weight_range = net.get_weights_range(layer_ind);
    p = frac_of_weight_range * layer_weight_range;
    net.Layers{layer_ind}.perturb_whole_layer(-p, p);
    result = net.verify_robustness_for_3dim_img(reachOptions, input = img);
    
    net_parallel = getArrayFromByteStream(netbs);
    net_parallel.Layers{layer_ind}.perturb_whole_layer(-p, p);
    result_parallel = net_parallel.verify_robustness_for_3dim_img(reachOptions_parallel, input = img);
    
    assert(result == result_parallel)
end
