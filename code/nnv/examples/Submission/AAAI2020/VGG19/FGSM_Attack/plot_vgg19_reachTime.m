
clc;
clear;


fprintf('\n\n=============================LOAD VGG19 ======================\n');

% Load the trained model 
net = vgg19();

fprintf('\n\n======================== PARSING VGG19 =======================\n');
nnvNet = CNN.parse(net, 'VGG19');


fprintf('\n\n=========CONSTRUCT INPUT SET (AN IMAGESTAR SET) =============\n');
load image_data.mat;
V(:,:,:,1) = double(ori_image);
V(:,:,:,2) = double(dif_image);

pred_lb = 0.5;
delta = 0.0000001;
pred_ub = pred_lb + delta;
C = [1;-1];   % pred_lb % <= alpha <= pred_ub percentage of FGSM attack
d = [pred_ub; -pred_lb];
IS = ImageStar(double(V), C, d, pred_lb, pred_ub);

fprintf('\n\n======= DO REACHABILITY ANLAYSIS WITH EXACT-STAR METHOD ======\n');

nnvNet.reach(IS, 'exact-star');

exactReachSet = nnvNet.reachSet{44};

fprintf('\n\n========REACHABILITY IS DONE IN %.5f SECONDS==========\n', nnvNet.totalReachTime);


RT = nnvNet.reachTime; 

conv_reachTime = RT(2) + RT(4) + RT(7) + RT(9) + RT(12) + RT(14) + RT(16) + RT(18) + RT(21) + RT(23) + RT(25) + RT(27) + RT(30) + RT(32) + RT(34) + RT(36);
fc_reachTime = RT(39) + RT(41) + RT(43);
mp_reachTime = RT(6) + RT(11) + RT(20) + RT(29) + RT(38);
rl_reachTime = RT(3) + RT(5) + RT(8) + RT(10) + RT(13) + RT(15) + RT(17) + RT(19) + RT(22) + RT(24) + RT(26) + RT(28) + RT(31) + RT(33) + RT(35) + RT(37) + RT(40) + RT(42);

figure;
c = categorical({'Convolutional 2D (16)','Fully Connected (3)', 'Max Pooling (5)', 'ReLU (18)'});
reachTime = [conv_reachTime fc_reachTime mp_reachTime rl_reachTime];
bar(c,reachTime);
text(1:length(reachTime),reachTime,num2str(reachTime'),'vert','bottom','horiz','center');




