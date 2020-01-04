% Load the trained model 
net = vgg16();
L1 = ImageInputLayer.parse(net.Layers(1));

% See details of the architecture 
net.Layers

% Read the image to classify 
I0 = imread('peppers.png');

% Adjust size of the image 
sz = net.Layers(1).InputSize; 
I = I0(1:sz(1),1:sz(2),1:sz(3));

n = size(I);
N = n(1)*n(2)*n(3);

I4 = reshape(I, [N,1]);
I5 = I4;
del = 200;

for i=1:N
    if I4(i) >= 255 - del
        I5(i) = 255 - I4(i);
    else
        I5(i) = 0;
    end
end

center = reshape(I4, [n(1), n(2), n(3)]); % center image matrix
basis_mat = reshape(I5, [n(1), n(2), n(3)]); % basis image matrix

V = double(cat(4, center, basis_mat));

tic;
image = ImageZono(V);
toc;

