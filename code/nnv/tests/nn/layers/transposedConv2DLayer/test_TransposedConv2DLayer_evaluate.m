% construct a Conv2DLayer object
% this convolution layer is from the link: 
% http://cs231n.github.io/convolutional-networks/

% filter 1 weight matrix with 3 channels
W(:,:, 1, 1) = [1 1 1; 1 1 -1; 0 -1 0]; % channel 1
W(:, :, 2, 1) = [0 1 1; 0 -1 0; 0 0 -1]; % channel 2

% filter 2 weight matrix with 3 channels
W(:, :, 1, 2) = [1 1 -1; -1 0 1; 1 -1 -1]; % channel 1
W(:, :, 2, 2) = [-1 1 1; 0 1 1; 0 -1 1]; % channel 2

% biases for 2 filters
b(:, :, 1) = 1; % filter 1
b(:, :, 2) = 0; % filter 2



L = TransposeConv2DLayer(W, b);
input = rand(4,4,2);
output = L.evaluate(input);
