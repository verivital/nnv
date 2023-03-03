path_out = [path_results(), filesep, 'vgg19', filesep];


fprintf('\n\n=============================LOAD VGG16 ======================\n');

clc;
clear;

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

pred_lb = 0.98;
delta = 0.0000002;
pred_ub = pred_lb + delta;
C = [1;-1];   % pred_lb % <= alpha <= pred_ub percentage of FGSM attack
d = [pred_ub; -pred_lb];
IS = ImageStar(double(V), C, d, pred_lb, pred_ub);

fprintf('\n\n========= PARSE VGG19 FOR REACHABILITY ANALYSIS ============\n');

nnvNet = CNN.parse(net, 'VGG19');

fprintf('\n\n======= DO REACHABILITY ANLAYSIS WITH EXACT-STAR METHOD ======\n');

nnvNet.reach(IS, 'exact-star');

exactReachSet = nnvNet.reachSet{44};

fprintf('\n\n========REACHABILITY IS DONE IN %.5f SECONDS==========\n', nnvNet.totalReachTime);



point1 = [1 1 946]; % bell pepper index
point2 = [1 1 950]; % strawberry index

S1 = [];
N = length(exactReachSet);
for i=1:N
    S1 = [S1 exactReachSet(i).project2D(point1, point2)];
end

fprintf('\n\n==================PLOT COUNTER EXAMPLE====================\n');


% load test image
load image_data.mat;

adv_image = double(ori_image) + pred_ub*double(dif_image);

label1 = classify(net, ori_image);
label2 = classify(net, adv_image);



figure; 
subplot(2,3,[1 2 3]);
Star.plotBoxes_2D(S1, 1, 2, 'blue');
xlabel('Bell Pepper (946)', 'FontSize', 13);
ylabel('Strawberry (950)', 'FontSize', 13);

subplot(2,3,4); 
imshow(ori_image);
text(10, 20, char(label1),'Color','black');
title('Original image');

subplot(2,3,5);
imshow(dif_image); 
title('Difference image');

subplot(2, 3, 6);
imshow((adv_image/255)); 
text(10, 20, char(label2),'Color','black');
xlabel('l = 98\%, $\delta = 2\times 10^{-7}$', 'Interpreter','latex', 'FontSize', 13);
title('Adversarial image');

path_out_vgg19 = [path_results(), filesep, 'vgg19', filesep];
saveas(gcf, [path_out_vgg19, 'figure10_vgg19.png'])

% there is a bug in the matlab scheduler
% adding this exit(0) avoids throwing an OOM error, which does not actually happen or cause any results to not be generated
% see here: https://www.mathworks.com/matlabcentral/answers/442711-script-fails-when-run-via-scheduler-matlab-management-cpp-671-find-assertion-failed
% note that with this exit though, this script must be called last, as it will kill the matlab session
exit(0);