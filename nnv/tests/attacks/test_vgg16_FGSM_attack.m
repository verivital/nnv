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
del1 = 0.1;
del2 = 0.2;
del3 = 0.3;

I31 = I2 + del1*r;
I32 = I2 + del2*r;
I33 = I2 + del3*r;

I31 = reshape(I31, [n(1), n(2), n(3)]);
I32 = reshape(I32, [n(1), n(2), n(3)]);
I33 = reshape(I33, [n(1), n(2), n(3)]);

figure; 
imshow(I31);
title("FGSM attack with del = 0.1");

figure; 
imshow(I32);
title("FGSM attack with del = 0.2");

figure; 
imshow(I33);
title("FGSM attack with del = 0.3");


% Classify the image using VGG-16 
label = classify(net, I);

% Show the image and the classification results 
figure; 
imshow(I) 
text(10, 20, char(label),'Color','white')


label1 = classify(net, I31);
label2 = classify(net, I32);
label3 = classify(net, I33);

% Show the image and the classification results 
figure; 
imshow(I31) 
text(10, 20, char(label1),'Color','white')

% Show the image and the classification results 
figure; 
imshow(I32) 
text(10, 20, char(label2),'Color','white')

% Show the image and the classification results 
figure; 
imshow(I33) 
text(10, 20, char(label3),'Color','white')

