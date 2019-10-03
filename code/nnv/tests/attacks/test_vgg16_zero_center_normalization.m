% Load the trained model 
net = vgg16();

% See details of the architecture 
net.Layers

% Read the image to classify 
I0 = imread('peppers.png');

% Adjust size of the image 
sz = net.Layers(1).InputSize; 
I = I0(1:sz(1),1:sz(2),1:sz(3));

n = size(I);
N = n(1)*n(2)*n(3);

% FGSM attack
I2 = reshape(I, [N,1]);

beta = 255;
r = randi(beta,N,1, 'uint8');
del1 = 0.2;
I31 = I2 + del1*r;

I31 = reshape(I31, [n(1), n(2), n(3)]);

figure; 
imshow(I31);
title("FGSM attack with del = 0.1");


% Classify the image using VGG-16 
label = classify(net, I);

% Show the image and the classification results 
figure; 
imshow(I) 
text(10, 20, char(label),'Color','white')


label1 = classify(net, I31);
% Show the image and the classification results 
figure; 
imshow(I31) 
text(10, 20, char(label1),'Color','white')

O31 = activations(net, I31, 1);
I32 = cast(reshape(I31,[N,1]), 'single');
O31 = reshape(O31, [N, 1]);


