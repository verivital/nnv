net = NN();

l = FullyConnectedLayer();
net.Layers{1} = l;
net.numLayers = 1;
l.InputSize = 2;
l.OutputSize = 2;
l.Weights = [1 0; 0 1];
l.Bias = [1; 0];
netbs = getByteStreamFromArray(net);
% net.Layers{1}

x = [1; 1];
y = net.evaluate(x);

I = ImageStar(x, x);
net.Layers{1}.perturb_whole_layer_given_fraction_of_weights_range(0.5);
O = net.reach(I);
O.getRanges();
lb = reshape(O.im_lb, [], 1);
ub = reshape(O.im_ub, [], 1);

% Sample the weight perturbed NN to compare bounds with ImageStar set
% computations. The results should be the same because the
% toy NN above has no non-linear layers.
[compare_bounds, mismatching_bounds] = WPutils.sample_weight_perturbed_nns(net, ...
                                input = x, ...
                                samples_per_pert = 2, ...
                                lb_out = lb, ...
                                ub_out = ub, ...
                                netbs = netbs, ...
                                robustness_result = 0, ...
                                tolerance_for_comparison_of_bounds = 1e-4);

assert(isempty(mismatching_bounds))