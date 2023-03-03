clc;
clear;

% Load mnist data set
load('test_images.mat');
% load cnn
load('TEST_NET.mat');

% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9
 
nnv_net = CNN.parse(net, 'TEST_CONVNET' );

IM = IM_data(:,:,1);
IM = reshape(IM, [28, 28, 1]);
imshow(IM);

[~,Y] = nnv_net.evaluate(IM);
Z9= activations(net,IM,9);
Z10 = activations(net, IM, 10);

dif1 = Y{8} - Z9;
max_dif1 = max(dif1, [], 'all');

