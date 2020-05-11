% construct a Conv2DLayer object
% this convolution layer is from the link: 
% http://cs231n.github.io/convolutional-networks/

% filter 1 weight matrix with 3 channels
W(:,:, 1, 1) = [1 1 1; 1 1 -1; 0 -1 0]; % channel 1
W(:, :, 1, 2) = [0 1 1; 0 -1 0; 0 0 -1]; % channel 2
W(:, :, 1, 3) = [1 -1 0; 1 1 0; -1 0 1]; % channel 3
% filter 2 weight matrix with 3 channels
W(:, :, 2, 1) = [1 1 -1; -1 0 1; 1 -1 -1]; % channel 1
W(:, :, 2, 2) = [-1 1 1; 0 1 1; 0 -1 1]; % channel 2
W(:, :, 2, 3) = [-1 1 -1; -1 0 -1; -1 -1 -1]; % channel 3

% biases for 2 filters
b(:, :, 1) = 1; % filter 1
b(:, :, 2) = 0; % filter 2


L0 = TransposeConv2DLayer(W, b);

LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image = ImageZono(LB, UB);

S = L0.reach(image.toImageStar);