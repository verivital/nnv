%to run this as a test, use results_layers_Conv2D=runtests('test_layers_Conv2D')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables




%___________________________________________________________________________________________________
%tests below originally taken from test_Conv2DLayer_compute_featureMap.m


%% test 1: Conv2DLayer feature map


I = [1 1 1 0 0; 0 1 1 1 0; 0 0 1 1 1; 0 0 1 1 0; 0 1 1 0 0]; % input
W = [1 0 1; 0 1 0; 1 0 1]; % filter

padding = [0 0 0 0];
stride = [1 1];
dilation = [1 1];


featureMap = Conv2DLayer.compute_featureMap(I, W, padding, stride, dilation);


checker=[I(1, 1)+I(1, 3)+I(2, 2)+I(3, 1)+I(3, 3), I(1, 2)+I(1, 4)+I(2, 3)+I(3, 2)+I(3, 4), I(1, 3)+I(1, 5)+I(2, 4)+I(3, 3)+I(3, 5);  I(2, 1)+I(2, 3)+I(3, 2)+I(4, 1)+I(4, 3), I(2, 2)+I(2, 4)+I(3, 3)+I(4, 2)+I(4, 4), I(2, 3)+I(2, 5)+I(3, 4)+I(4, 3)+I(4, 5); I(3, 1)+I(3, 3)+I(4, 2)+I(5, 1)+I(5, 3), I(3, 2)+I(3, 4)+I(4, 3)+I(5, 2)+I(5, 4), I(3, 3)+I(3, 5)+I(4, 4)+I(5, 3)+I(5, 5);];   

assert(isequal(checker, featureMap))




%___________________________________________________________________________________________________
%tests below originally taken from test_Conv2DLayer_constructor.m


%% test 2: Conv2DLayer constructor


% construct a Conv2DLayer object
% this convolution layer is from the link: 
% http://cs231n.github.io/convolutional-networks/

% filter 1 weight matrix with 3 channels
W(:,:, 1, 1) = [1 1 1; 1 1 -1; 0 -1 0]; % channel 1
W(:, :, 2, 1) = [0 1 1; 0 -1 0; 0 0 -1]; % channel 2
W(:, :, 3, 1) = [1 -1 0; 1 1 0; -1 0 1]; % channel 3
% filter 2 weight matrix with 3 channels
W(:, :, 1, 2) = [1 1 -1; -1 0 1; 1 -1 -1]; % channel 1
W(:, :, 2, 2) = [-1 1 1; 0 1 1; 0 -1 1]; % channel 2
W(:, :, 3, 2) = [-1 1 -1; -1 0 -1; -1 -1 -1]; % channel 3

% biases for 2 filters
b(:, :, 1) = 1; % filter 1
b(:, :, 2) = 0; % filter 2


L0 = Conv2DLayer(W, b);

L1 = Conv2DLayer('Test_Cov2DLayer', W, b, [1 1 1 1], [1 1], [1 1]);



%___________________________________________________________________________________________________
%tests below originally taken from test_Conv2DLayer_get_zero_padding_input.m


%% test 3: Conv2DLayer zero padding



% original input volume: color image with 3 channels
inputVol(:, :, 1) = [2 0 1 2 1; 1 0 2 2 2; 1 2 2 0 2; 1 2 0 0 1; 1 0 1 1 2]; % channel 1 input matrix
inputVol(:, :, 2) = [0 0 1 0 1; 0 0 2 1 1; 1 1 0 1 1; 1 1 0 2 2; 2 1 2 0 0]; % channel 2 input matrix
inputVol(:, :, 3) = [1 2 2 1 0; 2 0 0 2 0; 0 0 1 0 1; 1 2 0 2 0; 1 0 2 1 0]; % channel 3 input matrix


% construct input with padding operation
paddingSize = [1 1 1 1];
I = Conv2DLayer.get_zero_padding_input(inputVol, paddingSize);

assert(isequal(I(2:end-1, 2:end-1, :), inputVol))

corners=[I(1, 1, :), I(1, end, :); I(end, 1, :), I(end, end, :)];
top_bot=[I(1, 2:end-1, :); I(end, 2:end-1, :)];
left_right=[I(2:end-1, 1, :), I(2:end-1, end, :)];

assert(isempty(find(corners)))
assert(isempty(find(top_bot)))
assert(isempty(find(left_right)))



%___________________________________________________________________________________________________
%tests below originally taken from test_Conv2DLayer_reach.m


%% test 4: Conv2DLayer reach

% construct a Conv2DLayer object
% this convolution layer is from the link: 
% http://cs231n.github.io/convolutional-networks/

% filter 1 weight matrix with 3 channels
W(:,:, 1, 1) = [1 1 1; 1 1 -1; 0 -1 0]; % channel 1
W(:, :, 2, 1) = [0 1 1; 0 -1 0; 0 0 -1]; % channel 2
W(:, :, 3, 1) = [1 -1 0; 1 1 0; -1 0 1]; % channel 3
% filter 2 weight matrix with 3 channels
W(:, :, 1, 2) = [1 1 -1; -1 0 1; 1 -1 -1]; % channel 1
W(:, :, 2, 2) = [-1 1 1; 0 1 1; 0 -1 1]; % channel 2
W(:, :, 3, 2) = [-1 1 -1; -1 0 -1; -1 -1 -1]; % channel 3

% biases for 2 filters
b(:, :, 1) = 1; % filter 1
b(:, :, 2) = 0; % filter 2


L0 = Conv2DLayer(W, b);

LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image = ImageZono(LB, UB);

%Z = L0.reach(image); %THIS TRIGGERS AN ERROR, BECAUSE IT'S NOT AN IMAGE

S = L0.reach(image.toImageStar);



%___________________________________________________________________________________________________
%tests below originally taken from test_Conv2DLayer_reach_star_exact_single_input.m


%% test 5: Conv2DLayer reach star exact


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



















