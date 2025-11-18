
clear

nn_name = 'Small_ConvNet';
folder = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Submission', filesep, 'CAV2020_ImageStar', filesep, 'MNIST_NETS', filesep, 'Architecture', filesep];
mnist_model = load([folder nn_name '.mat']);
matlabnet = mnist_model.net;

% Create NNV model that has only the first two layers: ImageInputLayer and
% Conv2DLayer
original_net = matlab2nnv(matlabnet);
l1 = original_net.Layers{1};
l2 = original_net.Layers{2};
net = NN();
net.Layers = {l1, l2};
net.numLayers = 2;

netbs = getByteStreamFromArray(net);    % bytestream to copy this network
    % for weight perturbations in parfor loop iterations
        
% Define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approximate method
reachOptions.numCores = 1;
% reachOptions.device = 'gpu';
reachOptions.device = 'cpu';
% reachOptions.delete_old_sets = 1;
reachOptions.free_mem_frac_for_LPs = 0.1;
% reachOptions.dis_opt = 'display';
% reachOptions.disp_intersection_result = 1;
% reachOptions.debug = 1;

[images, labels] = load_images(database = "mnist", ...
    n = 10, ...
    matlabnet = matlabnet);
img = images{1};
V = img;
V(:,:,:,2) = zeros(size(img));
C = 0;
d = 0;
pred_lb = -1;
pred_ub = 1;
I = ImageStar(V, C, d, pred_lb, pred_ub);

conv_layer_ind = 2;
pert = 0.1;
net.Layers{conv_layer_ind}.add_pert([1 1 1 1], -pert, pert);

O = net.reach(I, reachOptions);
O.getRanges();
lb = reshape(O.im_lb, [], 1);
ub = reshape(O.im_ub, [], 1);

% Sample the weight perturbed NN to compare bounds with ImageStar set
% computations. The results should be the same because the
% toy NN above has no non-linear layers.
[compare_bounds, mismatching_bounds] = net.sample_weight_perturbed_nns(input = img, ...
                                samples_per_pert = 2, ...
                                lb_out = lb, ...
                                ub_out = ub, ...
                                netbs = netbs, ...
                                robustness_result = 0, ...
                                tolerance_for_comparison_of_bounds = 1e-4);

assert(isempty(mismatching_bounds))
