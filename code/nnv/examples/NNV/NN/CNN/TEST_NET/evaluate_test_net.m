clc;
clear;

% Load mnist data set
load('test_images.mat');
% load cnn
load('TEST_NET.mat');

% get images and its labels
N = 2000;  % get 100 images from test images

% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9
 
Y = zeros(N, 1);
% check if CNN work ok
correct_num = 0;
for i=1:N
    Y(i) = classify(net, IM_data(:,:,i));
    if Y(i) ~= IM_labels(i)
        fprintf('\nNetwork evaluations give wrong results for the %d^th image', i);
        fprintf('\nLabels = %d, classification result = %d', IM_labels(i), Y(i));
        fprintf('\nThe image is a digit %d', IM_labels(i) - 1);
        fprintf('\nThe classified digit is %d', Y(i));
    else
        correct_num = correct_num + 1;
        fprintf('\nNetwork evaluations give correct results for the %d^th image', i);
        fprintf('\nLabels = %d, classification result = %d', IM_labels(i), Y(i));
        fprintf('\nThe image is a digit %d', Y(i) - 1);
    end
end

fprintf('\nClassification correctness: %.5f', correct_num/N);