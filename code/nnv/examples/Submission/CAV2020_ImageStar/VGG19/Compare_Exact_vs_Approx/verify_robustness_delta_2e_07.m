path_out = [path_results(), filesep, 'vgg19', filesep];
if ~isfolder(path_out)
    mkdir(path_out);
end

fprintf('\n\n=============================LOAD VGG19 ======================\n');

% Load the trained model 
if is_codeocean()
    load('/data/vgg19_cache.mat');
    net = net_vgg19;
else
    net = vgg19();
end

fprintf('\n\n======================== PARSING VGG19 =======================\n');
nnvNet = CNN.parse(net, 'VGG19');


fprintf('\n\n=========CONSTRUCT INPUT SET (AN IMAGESTAR SET) =============\n');
load image_data.mat;
V(:,:,:,1) = double(ori_image);
V(:,:,:,2) = double(dif_image);

l = [0.5 0.8 0.95 0.97 0.98 0.98999];
delta = 0.0000002;
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


fprintf('\n\n===================================VERIFICATION RESULTS===================================\n\n');

fprintf('\n l                    delta               Exact Analysis          Approximate Analysis');
fprintf('\n                                         Robust        VT           Robust         VT');
for i=1:n
    fprintf('\n %.6f          %.8f              %d        %.5f           %d         %.5f', l(i), delta, robust_exact(i), VT_exact(i), robust_approx(i), VT_approx(i));
end


save([path_out, 'table5_vgg19_verificationResult_2e_07.mat'], 'robust_exact', 'VT_exact', 'robust_approx', 'VT_approx');






