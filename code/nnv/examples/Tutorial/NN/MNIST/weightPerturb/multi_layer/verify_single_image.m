%% Robustness verification of a NN (L infinity adversarial attack)
%  if f(x) = y, then forall x' in X s.t. ||x - x'||_{\infty} <= eps,
%  then f(x') = y = f(x)

clear

% Load network
folder = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, 'NN', filesep, 'MNIST', filesep, 'weightPerturb', filesep];
mnist_model = load([folder 'mnist_model_fc.mat']);
classes = string(1:10);
matlabnet = mnist_model.net;

% Create NNV model
net = matlab2nnv(matlabnet);
net.OutputSize = length(classes);
netbs = getByteStreamFromArray(net);    % bytestream to copy this network
    % for weight perturbations in parfor loop iterations

[images, labels] = load_images(database = "mnist", ...
    n = 20, ...
    matlabnet = matlabnet);

target = single(labels{10}); % label = 4 (index 5 for our network)
img = single(images{10}); % change precision
correct_output = net.evaluate(img);

% % Visualize image;
% figure(1)
% imshow(img);

l1 = 7;
l2 = 9;
l3 = 13;

f = 0.0001;
p = f*net.get_weights_range(l1);
net.Layers{l1}.perturb_whole_layer(-p, p);
net.Layers{l2}.perturb_whole_layer(-p, p);

reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approximate method
reachOptions.numCores = 1;
reachOptions.device = 'cpu';
reachOptions.free_mem_frac_for_LPs = 0.1;
% reachOptions.layer_specific_numCores = dictionary( ...
%     ["start", "ReluLayer", "MaxPooling2DLayer", "end"], ...
%     [  64   ,     64     ,         64         ,  64  ]);
% parpool(parcluster('local').NumWorkers);

% Verification
t = tic;

result = net.verify_robustness_for_3dim_img(reachOptions, input = img); % target not provided because
    % the images to be verified are selected for verification
    % in load_images only if the NN classifies them correctly
    % in the absence of perturbations.

if result == 1
    disp("Neural network is verified to be robust!")
% elseif result == 0
%     disp("Neural network is NOT robust!")
else
    disp("Unknown result")
end

toc(t);


%% Let's visualize the ranges for every possible output

% Get output reachable set
R = net.reachSet{end};

check_bounds_match = 0;
if length(R) == 1
    % Get (overapproximate) ranges for each output index
    [imagestar_lb, imagestar_ub] = R.getRanges;
    imagestar_lb = squeeze(imagestar_lb);
    imagestar_ub = squeeze(imagestar_ub);
    
    % Get middle point for each output and range sizes
    mid_range = (imagestar_lb + imagestar_ub)/2;
    range_size = imagestar_ub - mid_range;
    
    % Label for x-axis
    x = [0 1 2 3 4 5 6 7 8 9];
    
    % % Visualize set ranges and evaluation points
    % figure(2);
    % hold off;
    % errorbar(x, mid_range, range_size, '.');
    % hold on;
    % xlim([-0.5 9.5]);
    % scatter(x,correct_output, 'x', 'MarkerEdgeColor', 'r');
else
    check_bounds_match = 0;
    
    figure(5);
    tiledlayout(1, length(R));
    imagestar_lb = inf(size(correct_output));
    imagestar_ub = -imagestar_lb;
    for IS_no=1:length(R)
        % Get (overapproximate) ranges for each output index
        nexttile
        [lb, ub] = R(IS_no).getRanges;
        lb = squeeze(lb);
        ub = squeeze(ub);
        
        imagestar_lb = min(imagestar_lb, lb);
        imagestar_ub = max(imagestar_ub, ub);
        
        % Get middle point for each output and range sizes
        mid_range = (lb + ub)/2;
        range_size = ub - mid_range;
        
        % Label for x-axis
        x = [0 1 2 3 4 5 6 7 8 9];
        
        % Visualize set ranges and evaluation points
        hold off;
        errorbar(x, mid_range, range_size, '.');
        hold on;
        xlim([-0.5 9.5]);
        scatter(x,correct_output, 'x', 'MarkerEdgeColor', 'r');
    end
end

%%
% Towards Certificated ...

xN = net.input_sets{l1}.V(:,:,:,1);
xN = xN(:);
yN = net.Layers{l1}.Weights*xN + net.Layers{l1}.Bias;

% perturbing whole layer by [-p,+p]
pert_vec = p*ones(1, length(xN));

% % additional change if only perturbing half layer by [-p,+p]
% pert_vec(end/2 + 1:end) = zeros(1, length(pert_vec(end/2 + 1:end)));

% get bounds of the layer before l2 because l2 will be perturbed as well
[pre_activation_lb_before_l2, pre_activation_ub_before_l2] = compute_bounds_weng_2020(net, pert_vec, inf, l1, l2 - 2);
lb_before_l2 = max(0, pre_activation_lb_before_l2);
ub_before_l2 = max(0, pre_activation_ub_before_l2);


[towards_lb, towards_ub] = propagate_bounds_weng_2020(net, p, l2, l3, lb_before_l2, ub_before_l2);
towards_mid = (towards_lb + towards_ub)/2;
towards_range = towards_ub - towards_mid;

% if towards_lb < imagestar_lb + abs(imagestar_lb)/1e4 & towards_ub > imagestar_ub - abs(imagestar_ub)/1e4
% else
%     error('Potential problem: bounds seem inconsistent!')
% end

%% Visualize set ranges and evaluation points
figure(2);
tiledlayout(1,2, 'TileSpacing', 'compact');
ax1 = nexttile;
hold off;
% imagestar
errorbar(x, mid_range, range_size, '.');
hold on;
xlim([-0.5 9.5]);
scatter(x,correct_output, 'x', 'MarkerEdgeColor', 'r');
xticks(0:9)
xlabel('MNIST Classes')
ylabel('Over-approximated Ranges of Output Layer')
title({'Modelstar', '\rm({\bfVerified} Robustness)'})

ax2 = nexttile;
hold off;
% towards ...
errorbar(x, towards_mid, towards_range, '.');
hold on;
xlim([-0.5 9.5]);
scatter(x,correct_output, 'x', 'MarkerEdgeColor', 'r');
xticks(0:9)
xlabel('MNIST Classes')
yticklabels([]);
title({'Certificated-Robust', '\rm(Did {\bfNot Verify} Robustness)'})
linkaxes([ax1 ax2], 'y')

lg = legend({'Output Reachable Set Ranges', 'Unperturbed Output'}, 'Orientation', 'horizontal');
lg.Layout.Tile = 'North';
% lg.Box = 'off';

fig = gcf;
fig.Position = [398 116 512 413];

%% function for "Towards..."
function [lb, ub] = compute_bounds(net, pert_vec, N, K)
    % N is the perturbed layer
    % K is the layer whose output's bounds are returned
    
    layers = [];
    for k = N:length(net.Layers)
        if isa(net.Layers{k}, 'FullyConnectedLayer')
            layers = [layers k];
        end
        if k == K
            break
        end
    end
    
    xN = net.input_sets{N}.V(:,:,:,1);
    xN = xN(:);
    yN = net.Layers{N}.Weights*xN + net.Layers{N}.Bias;
    pert = pert_vec*abs(xN);
    lb = yN - pert;
    ub = yN + pert;
    for k = 2:length(layers)
        bs = [lb ub];
        bs(bs < 0) = 0;
        W = net.Layers{layers(k)}.Weights;
        b = net.Layers{layers(k)}.Bias;
        lb = zeros(size(W, 1), 1);
        ub = lb;
        for m = 1:size(W, 1)
            rows = W(m,:).*bs';
            lb(m) = sum(min(rows)) + b(m);
            ub(m) = sum(max(rows)) + b(m);
        end
    end
end


function [lb, ub] = compute_bounds_only_used_slope(net, pert_vec, N, K)
    % N is the perturbed layer
    % K is the layer whose output's bounds are returned
    
    layers = [];
    for k = N:length(net.Layers)
        if isa(net.Layers{k}, 'FullyConnectedLayer')
            layers = [layers k];
        end
        if k == K
            break
        end
    end
    
    xN = net.input_sets{N}.V(:,:,:,1);
    xN = xN(:);
    yN = net.Layers{N}.Weights*xN + net.Layers{N}.Bias;
    pert = pert_vec*abs(xN);
    lb = yN - pert;
    ub = yN + pert;
    for k = 2:length(layers)
        linear_lb = ub./(ub - lb).*lb;
        lb(lb < 0) = linear_lb(lb < 0);
        lb(ub < 0) = 0;
        ub(ub < 0) = 0;
        W = net.Layers{layers(k)}.Weights;
        b = net.Layers{layers(k)}.Bias;
        lbk = zeros(size(W, 1), 1);
        ubk = lbk;
        for m = 1:size(W, 1)
            rows = W(m,:).*[lb ub]';
            lbk(m) = sum(min(rows)) + b(m);
            ubk(m) = sum(max(rows)) + b(m);
        end
        lb = lbk;
        ub = ubk;
    end
end

