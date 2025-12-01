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

% transform into NNV
net = matlab2nnv(matlab_net);
net.Name = 'collins_rul_cnn_small_window_20';
% net.OutputSize = 1;
netbs = getByteStreamFromArray(net);    % bytestream to copy this network
    % for weight perturbations in parfor loop iterations

% Define reachability parameters
reachOptions.reachMethod = 'approx-star';
reachOptions.device = 'cpu';

% Perturb weights
conv2d_layer_indices = net.get_layer_indices(["Conv2DLayer"]);
perturbed_layer_index = conv2d_layer_indices(1);

f = 0.05;	% fraction of layer weight range to use as perturbation
% At f = 0.02, the property is satisfied
% At f = 0.05, the property is not satisfied

% copy neural network from blank (unmodified)
net = getArrayFromByteStream(netbs);

net.Layers{perturbed_layer_index}.perturb_whole_layer_given_fraction_of_weights_range(f);

% Verify network
t = tic;
result = net.verify_vnnlib(vnnlib_path, reachOptions, needReshape);
time = toc(t);
if result == 1
    disp("Property satisfied")
else
    disp("Property NOT satisfied")
end
