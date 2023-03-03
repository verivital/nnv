%% parse network into NNV
modelfile = 'model_MNIST_CNN_Small.json';
weightfile = 'model_MNIST_CNN_Small.h5';

% modelfile = 'model_MNIST_CNN_Medium.json';
% weightfile = 'model_MNIST_CNN_Medium.h5';

% modelfile = 'model_MNIST_CNN_Large.json';
% weightfile = 'model_MNIST_CNN_Large.h5';

net = importKerasNetwork(modelfile, 'WeightFile', weightfile, 'OutputLayerType','classification');
nnvNet = CNN.parse(net); % construct an nnvNet object

%% Load an image and create attack
%load test_images.mat; 
%im = im_data(:,:,1);
%im = im/255;

load digit_7.mat; % the first image in JiaMeng data set
im = digit_7/255;
im = im';


% attack all pixels independently by some bounded disturbance d 
% d = 0.001; % 63 seconds
% d = 0.002; % 100 seconds
% d = 0.003; % 160 seconds
% d = 0.004; % 214 seconds
d = 0.01; 
attack_LB = -d*ones(28, 28);
attack_UB = d*ones(28,28);

%% computing reachable set 

IS = ImageStar(im, attack_LB, attack_UB); % construct an ImageStar input set
% t = tic;
% OS1 = nnvNet.reach(IS, 'approx-star'); % perfrom reachability analysis using approx-star method
% OS1.estimateRanges; % estimate ranges of all outputs
% rt1 = toc(t);

t = tic;
relaxFactor = 0.9;
OS2 = nnvNet.reach(IS, 'approx-star', 1, relaxFactor); % perfrom reachability analysis using approx-star method
OS2.estimateRanges; % estimate ranges of all outputs
rt2 = toc(t);

%% get output ranges
% output ranges with ImageStar method
% lb1 = reshape(OS1.im_lb, [10 1]); % lower bound vector of the output
% ub1 = reshape(OS1.im_ub, [10 1]); % upper bound vector of the output

lb2 = reshape(OS2.im_lb, [10 1]); % lower bound vector of the output
ub2 = reshape(OS2.im_ub, [10 1]); % upper bound vector of the output

%[lb_8, ub_8] = OS1.getRange(1,1,8); % digit 7

% lb2 = reshape(lb2, [10 1]);
% ub2 = reshape(ub2, [10 1]);

% %% plot output ranges
% 
% im_center1 = (lb1 + ub1)/2;
% err1 = (ub1 - lb1)/2;
% x1 = 0:1:9;
% y1 = im_center1;
% 
% figure;
% 
% e = errorbar(x1,y1,err1);
% e.LineStyle = 'none';
% e.LineWidth = 1;
% e.Color = 'red';
% xlabel('Output', 'FontSize', 11);
% ylabel('Ranges', 'FontSize', 11);
% xlim([0 9]);
% title('ImageStar', 'FontSize', 11);
% xticks(x1);
% xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'});
% set(gca, 'FontSize', 10);

