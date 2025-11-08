%% Script to verify some properties for Collins RUL CNN

clear

rul_path = ['.'];
needReshape = 2;

nn_path = './NN_rul_small_window_20.onnx';
vnnlib_path = './monotonicity_CI_shift10_w20.vnnlib';

reachOptions = struct;
reachOptions.single_average_input = 1;

% Select network to verify
matlab_net = importNetworkFromONNX(nn_path);

fprintf("Network: %s\n", net_name);

% transform into NNV
net = matlab2nnv(matlab_net);
net.Name = 'collins_rul_cnn_small_window_20';
% net.OutputSize = 1;
netbs = getByteStreamFromArray(net);    % bytestream to copy this network
    % for weight perturbations in parfor loop iterations

% Define reachability parameters
reachOptions.reachMethod = 'approx-star';
reachOptions.device = 'cpu';
reachOptions.delete_old_sets = 1;
reachOptions.free_mem_frac_for_LPs = 0.1;
% reachOptions.dis_opt = 'display';
% reachOptions.disp_intersection_result = 1;

% designation of number of workers to use for parallelization; skip the "layer_specific_numCores" field to avoid parallelization within layers altogether
reachOptions.layer_specific_numCores = dictionary( ...
    ["start", "ReluLayer", "MaxPooling2DLayer", "end"], ...
    [  64   ,     64     ,         64         ,  64  ]);	% fix the number of workers to use; "end" refers to verification stage

% reachOptions.layer_specific_numCores = dictionary( ...
%     ["start", "ReluLayer", "MaxPooling2DLayer", "end"], ...
%     [  64   ,     -1     ,         -1         ,  -1  ]);	% -1 means, estimate the appropriate number of workers using regression-based polynomial for memory usage estimation for star set LPs

% Perturb weights
conv2d_layer_indices = net.get_layer_indices(["Conv2DLayer"]);
perturbed_layer_index = conv2d_layer_indices(1);

f = 1e-3;	% fraction of layer weight range to use as perturbation

% copy neural network from blank (unmodified)
net = getArrayFromByteStream(netbs);

net.Layers{l_index_in_net}.perturb_whole_layer_given_fraction_of_weights_range(f);

% Verify network
t = tic;
modelstar_result = net.verify_vnnlib(vnnlib_path, reachOptions, needReshape);
modelstar_time = toc(t);
if modelstar_result == 2
    modelstar_result = 0;
end

fprintf("m%d(%.3fs) c%d(%.3fs) ", modelstar_result, modelstar_time, cert_robust_result, cert_robust_time);



