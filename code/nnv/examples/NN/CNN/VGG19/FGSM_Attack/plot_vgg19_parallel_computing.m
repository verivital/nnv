
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
delta = 0.00000018;
pred_ub = pred_lb + delta;
C = [1;-1];   % pred_lb % <= alpha <= pred_ub percentage of FGSM attack
d = [pred_ub; -pred_lb];
IS = ImageStar(double(V), C, d, pred_lb, pred_ub);


numCores = [1 2 3 4];
n = length(numCores);
reachTime = zeros(n, 1);
for i=1:n
    nnvNet.reach(IS, 'exact-star', numCores(i));
    reachTime(i) = nnvNet.totalReachTime;
end

figure; 
plot(numCores, reachTime, '--o');
ax = gca;
ax.FontSize = 13; 
xlabel('Number of Cores','FontSize',13);
ylabel('Reachability Time', 'FontSize', 13);



