%to run this as a test, use results_layers_TransposedConv2DLayer=runtests('test_layers_TransposedConv2DLayer')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables




%___________________________________________________________________________________________________
%tests below originally taken from test_TransposedConv2DLayer_evaluate.m


%% test 1: TransposedConv2DLayer evaluate

% construct a Conv2DLayer object
% this convolution layer is from the link: 
% http://cs231n.github.io/convolutional-networks/

% filter 1 weight matrix with 3 channels
W(:,:, 1, 1) = [1 1 1; 1 1 -1; 0 -1 0]; % channel 1
W(:, :, 2, 1) = [0 1 1; 0 -1 0; 0 0 -1]; % channel 2
%W(:, :, 2, 1) = [0 0 0; 0 0 0; 0 0 0]; % channel 2

% filter 2 weight matrix with 3 channels
W(:, :, 1, 2) = [1 1 -1; -1 0 1; 1 -1 -1]; % channel 1
W(:, :, 2, 2) = [-1 1 1; 0 1 1; 0 -1 1]; % channel 2
%W(:, :, 1, 2) = [0 0 0; 0 100 0; 0 0 0]; % channel 1
%W(:, :, 2, 2) = [0 0 0; 0 0 0; 0 0 0]; % channel 2



% biases for 2 filters
b(:, :, 1) = 1; % filter 1
b(:, :, 2) = 0; % filter 2



L = TransposedConv2DLayer(W, b);
input = rand(4,4,2);
%input(:, :, 1)=[1, 2, 3, 4; 5, 6, 7, 8; 1, 2, 3, 4; 1, 1, 1, 1];
%input(:, :, 1)=1000*[1, 2, 3, 4; 5, 6, 7, 8; 1, 2, 3, 4; 1, 1, 1, 1];
%input(:, :, 2)=[1, 2, 3, 4; 5, 6, 7, 8; 1, 2, 3, 4; 1, 1, 1, 1];
output = L.evaluate(input);

i2=zeros(8, 8, 2);
i2(3:6, 3:6, :)=input;

w_true=flip(flip(W, 1), 2);
conclusion=zeros(6, 6, 8);
for i=1:6
    for j=1:6
        for k1=1:2
            for k2=1:2
                for k3=1:2
                    conclusion(i, j, (k1-1)*4+(k2-1)*2+k3)=sum(i2(i:i+2, j:j+2, k1).*w_true(:, :, k2, k3), 'all');
                end
            end
        end
    end
end
o1=conclusion(:, :, 1)+conclusion(:, :, 6)+b(:, :, 1);
%this means what
%1-> i2(:, :, 1).*w_true(:, :, 1, 1)
%6-> i2(:, :, 2).*w_true(:, :, 1, 2)

o2=conclusion(:, :, 3)+conclusion(:, :, 8)+b(:, :, 2);
%this means what
%3-> i2(:, :, 1).*w_true(:, :, 2, 1)
%8-> i2(:, :, 2).*w_true(:, :, 2, 2)


%in particular, this means that each of these inputs changes both outputs.

assert(isequal(o1, output(:, :, 1)))
assert(isequal(o2, output(:, :, 2)))
%THESE CURRENTLY FAIL BECAUSE OF ROUNDING ERRORS.

%___________________________________________________________________________________________________
%tests below originally taken from test_TransposedConv2DLayer_reach.m


%% test 2: TransposedConv2DLayer reach

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


L0 = TransposedConv2DLayer(W, b);

LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image = ImageZono(LB, UB);

S = L0.reach(image.toImageStar);
