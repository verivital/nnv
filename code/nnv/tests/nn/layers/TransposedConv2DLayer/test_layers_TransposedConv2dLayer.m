%% test 1: TransposedConv2DLayer evaluate

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

L = TransposedConv2DLayer(W, b);
input = rand(4,4,2);
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
o2=conclusion(:, :, 3)+conclusion(:, :, 8)+b(:, :, 2);

assert(all(abs(o1 - output(:, :, 1)) <= 1e-12, 'all')); % 1e-12 is the tolerance for rounding errors
assert(all(abs(o2 - output(:, :, 2)) <= 1e-12, 'all'));


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
