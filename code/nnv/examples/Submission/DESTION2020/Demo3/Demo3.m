%% Using NNV to analyze the robustness of CNNs (VGG16 - IMAGENET)

%% Load the CNN (VGG16) and data
% Load the network
net = vgg16();
% Parse the network
nnvNet = CNN.parse(net, 'VGG16');

% Load image
load('image_data_vgg.mat');
V(:,:,:,1) = double(ori_image);
V(:,:,:,2) = double(dif_image);

%% Setup analysis and initialize variables
% l = [0.5 0.8 0.95 0.97 0.98 0.98999];
l = 0.5;
% l = [0.5 0.8 0.95];
delta = 0.0000001;
% delta = 0.0000002;
% This combinations shows the difference between exact and approach
% In both cases we show the NN is robust.
% Exact -> 801.49318 s   |     Approx -> 58.93040 s
% l = 0.9;
% delta = 0.00001;
% In this case my computer runs out of memory
% l = 0.98999;
% delta = 0.005;
n = length(l);

pred_lb = zeros(n, 1);
pred_ub = zeros(n, 1);
robust_exact = zeros(n, 1);
robust_approx = zeros(n, 1);
VT_exact = zeros(n, 1);
VT_approx = zeros(n, 1);

correct_id = 946; % id for bell pepper
for i=1:n
    pred_lb(i) = l(i);
    pred_ub(i) = l(i) + delta;
    
    C = [1;-1];   % pred_lb % <= alpha <= pred_ub percentage of FGSM attack
    d = [pred_ub(i); -pred_lb(i)]; 
    IS = ImageStar(V, C, d, pred_lb(i), pred_ub(i));
    
    fprintf('\n\n=================== VERIFYING ROBUSTNESS FOR l = %.5f, delta = %.8f ====================\n', l(i), delta);
    t = tic;
    [robust_exact(i), ~] = nnvNet.verifyRobustness(IS, correct_id, 'exact-star');
    VT_exact(i) = toc(t);
    t = tic; 
    [robust_approx(i), ~] = nnvNet.verifyRobustness(IS, correct_id, 'approx-star');
    VT_approx(i) = toc(t);
end

%% Show results
n = length(l);
fprintf('\n\n===================================VERIFICATION RESULTS===================================\n\n');

fprintf('\n l                          delta                        Exact Analysis          Approximate Analysis');
fprintf('\n                                                            Robust        VT             Robust         VT');
for i=1:n
    fprintf('\n %.6f          %.8f              %d        %.5f           %d         %.5f', l(i), delta, robust_exact(i), VT_exact(i), robust_approx(i), VT_approx(i));
end

%% Save results
% save verificationResult_1e_07.mat l delta robust_exact VT_exact robust_approx VT_approx;
% save verificationResult_L095_d1e_03.mat l delta robust_exact VT_exact robust_approx VT_approx;

