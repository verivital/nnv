%% parse network into NNV
modelfile = 'models/model_MNIST_CNN_Small.json';
weightfile = 'models/model_MNIST_CNN_Small.h5';

% Load keras model into MATLAB
net = importKerasNetwork(modelfile, 'WeightFile', weightfile, 'OutputLayerType','classification');
% convert MATLAB model into NNV
nnvNet = matlab2nnv(net); % construct an NN object

%% Load an image and create attack

load digit_7.mat; % the first image in JiaMeng data set
im = digit_7/255; % normalize data
im = im';

% attack all pixels independently by some bounded disturbance d 
d = 0.01; 
attack_LB = -d*ones(28, 28);
attack_UB = d*ones(28,28);

%% computing reachable set 

% Create input set as an ImageStar
IS = ImageStar(im, attack_LB, attack_UB); % construct an ImageStar input set

% Define reachability options
reachOptions = struct;
reachOptions.reachMethod = 'relax-star-area';
reachOptions.relaxFactor = 0.9;

% Compute reachability
t = tic;
OS2 = nnvNet.reach(IS, reachOptions); % perfrom reachability analysis using relax-star method
toc(t);

%% get output ranges
% output ranges with ImageStar method

[lb, ub] = OS2.getRanges(); % max should be index 8 (corresponds to label7)

% plot output ranges

im_center = (lb + ub)/2;
err1 = (ub - lb)/2;
x1 = 0:1:9;
y1 = im;

figure;
e = errorbar(x1,y1,err1);
e.LineStyle = 'none';
e.LineWidth = 1;
e.Color = 'red';
% xlabel('Output', 'FontSize', 11);
% ylabel('Ranges', 'FontSize', 11);
xlim([0 9]);
% title('ImageStar', 'FontSize', 11);
xticks(x1);
xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'});
% set(gca, 'FontSize', 10);

