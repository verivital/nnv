% construct a Conv2DLayer object
% this convolution layer is from the link: 
% http://cs231n.github.io/convolutional-networks/

L = Conv2DLayer(3, 2);
L.set_name('Test_Convolutional_Layer');

% filter has a size of 3 x 3
% filter 1 weight matrix with 3 channels
W(:,:, 1, 1) = [0 1 0; -1 1 -1; 1 1 1]; % channel 1
W(:, :, 2, 1) = [-1 -1 -1; 1 1 0; 0 0 1]; % channel 2
W(:, :, 3, 1) = [1 1 0; 1 1 1; -1 0 1]; % channel 3
% filter 2 weight matrix with 3 channels
W(:, :, 1, 2) = [1 -1 0; -1 1 1; 1 -1 1]; % channel 1
W(:, :, 2, 2) = [1 0 0; 0 1 -1; 0 -1 0]; % channel 2
W(:, :, 3, 2) = [-1 0 0; -1 -1 -1; 0 1 -1]; % channel 3

% biases for 2 filters
b(:, :, 1) = 1; % filter 1
b(:, :, 2) = 0; % filter 2


L.set_weights_biases(W, b);


% the convolutional layer has 2 filters of size 3x3
padding = 1;
stride = 2;
dilation = 1;

L.set_stride(stride);
L.set_padding(padding);
L.set_dilation(dilation);

% original input volume: color image with 3 channels
inputVol(:, :, 1) = [0 1 0 1 2; 2 2 0 1 2; 2 0 1 0 1; 0 2 2 2 1; 2 2 1 0 2]; % channel 1 input matrix
inputVol(:, :, 2) = [2 1 0 2 0; 1 0 2 0 2; 1 0 0 2 1; 0 2 0 2 1; 1 1 2 1 2]; % channel 2 input matrix
inputVol(:, :, 3) = [1 1 1 0 0; 2 0 1 0 0; 2 1 2 0 0; 2 0 0 0 1; 0 1 2 0 0]; % channel 3 input matrix


y = L.evaluate(inputVol);






