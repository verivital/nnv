% construct a Conv2DLayer object
% this convolution layer is from the link: 
% http://cs231n.github.io/convolutional-networks/

% filter has a size of 3 x 3
% filter 1 weight matrix with 3 channels

W(:,:, 1, 1) = [-1 0 0; 1 1 0; -1 0 -1]; % channel 1
W(:, :, 2, 1) = [-1 0 -1; 1 -1 1; 1 1 0]; % channel 2
W(:, :, 3, 1) = [-1 1 1; -1 -1 -1; 0 0 -1]; % channel 3
% filter 2 weight matrix with 3 channels
W(:, :, 1, 2) = [-1 1 0; 0 0 -1; 1 0 0]; % channel 1
W(:, :, 2, 2) = [1 0 0; 0 0 -1; 0 0 1]; % channel 2
W(:, :, 3, 2) = [1 -1 -1; 1 1 0; -1 0 1]; % channel 3

% biases for 2 filters
b(:, :, 1) = 1; % filter 1
b(:, :, 2) = 0; % filter 2

L = Conv2DLayer(W, b);
L.set_weights_biases(W, b);


% the convolutional layer has 2 filters of size 3x3
padding = 1;
stride = 2;
dilation = 1;

L.set_stride(stride);
L.set_padding(padding);
L.set_dilation(dilation);

% original input volume: color image with 3 channels
IM(:, :, 1) = [0 0 2 0 0; 1 2 0 2 0; 0 0 2 2 0; 0 2 2 2 2; 2 2 2 1 1]; % channel 1 input matrix
IM(:, :, 2) = [1 2 2 1 2; 2 1 2 0 2; 2 2 2 0 1; 1 1 1 0 0; 1 0 2 2 1]; % channel 2 input matrix
IM(:, :, 3) = [0 0 2 2 1; 0 2 1 1 2; 0 2 0 0 1; 0 2 1 0 1; 1 2 1 0 0]; % channel 3 input matrix


LB(:,:,1) = [-0.1 -0.2 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.2 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0;0 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0;0 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

input = ImageStar(IM, LB, UB);

Y = L.reach_star_single_input(input);






