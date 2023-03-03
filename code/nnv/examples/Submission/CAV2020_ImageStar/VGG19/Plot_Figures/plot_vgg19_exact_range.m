fprintf('\n\n=============================LOAD VGG19 ======================\n');

% Load the trained model 
if is_codeocean()
    load('/data/vgg19_cache.mat');
    net = net_vgg19;
else
    net = vgg19();
end

fprintf('\n\n=========CONSTRUCT INPUT SET (AN IMAGESTAR SET) =============\n');
load image_data.mat;
V(:,:,:,1) = double(ori_image);
V(:,:,:,2) = double(dif_image);

pred_lb = 0.95;
delta = 0.0000002;
pred_ub = pred_lb + delta;
C = [1;-1];   % pred_lb % <= alpha <= pred_ub percentage of FGSM attack
d = [pred_ub; -pred_lb];
IS = ImageStar(double(V), C, d, pred_lb, pred_ub);

fprintf('\n\n========= PARSE VGG19 FOR REACHABILITY ANALYSIS ============\n');

net = CNN.parse(net, 'VGG19');

fprintf('\n\n======= DO REACHABILITY ANLAYSIS WITH EXACT-STAR METHOD ======\n');

net.reach(IS, 'exact-star');

exactReachSet = net.reachSet{44};

fprintf('\n\n========REACHABILITY IS DONE IN %.5f SECONDS==========\n', net.totalReachTime);



point1 = [1 1 946]; % bell pepper index
point2 = [1 1 950]; % strawberry index

S1 = [];
N = length(exactReachSet);
for i=1:N
    S1 = [S1 exactReachSet(i).project2D(point1, point2)];
end

point1 = [1 1 946]; % bell pepper index
point2 = [1 1 807]; % sock index

S2 = [];
for i=1:N
    S2 = [S2 exactReachSet(i).project2D(point1, point2)];
end


im_lb = [];
im_ub = [];
for i=1:N
    [lb, ub] = exactReachSet(i).getRanges;
    lb = reshape(lb, [1000, 1]);
    ub = reshape(ub, [1000, 1]);
    im_lb = [im_lb lb];
    im_ub = [im_ub ub];
end

im_lb = min(im_lb, [], 2);
im_ub = max(im_ub, [], 2);


fprintf('\n\n===========================PLOT OUTPUT RANGES===============================\n')

im_center = (im_lb + im_ub)/2;
err = (im_ub - im_ub)/2;

x = 1:1:1000;
y = im_center;

figure;

subplot(2,2,[1 2]);
e = errorbar(x,y,err);
e.LineStyle = 'none';
e.Color = 'red';
xlabel('Output Category ID', 'FontSize', 13);
ylabel('Range', 'FontSize', 13);


subplot(2,2,3)
Star.plotBoxes_2D(S2, 1, 2, 'blue');
xlabel('Bell Pepper (946)', 'FontSize', 13);
ylabel('Sock (807)', 'FontSize', 13);


subplot(2,2,4)
Star.plotBoxes_2D(S1, 1, 2, 'blue');
xlabel('Bell Pepper (946)', 'FontSize', 13);
ylabel('Strawberry (950)', 'FontSize', 13);
