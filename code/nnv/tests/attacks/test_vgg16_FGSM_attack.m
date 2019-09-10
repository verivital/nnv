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
I21 = double(I2);

beta = 255;
%r = beta*rand(N,1);
del1 = 0.1;
del2 = 0.1593543;
del3 = 0.5;

I31 = I21 + del1*r;
I32 = I21 + del2*r;
I33 = I21 + del3*r;

I31 = reshape(I31, [n(1), n(2), n(3)]);
I32 = reshape(I32, [n(1), n(2), n(3)]);
I33 = reshape(I33, [n(1), n(2), n(3)]);


% Show the image and the classification results 
figure; 
subplot(1,4,1);
imshow(I);
title('Original image');



label1 = classify(net, I31);
label2 = classify(net, I32);
label3 = classify(net, I33);

% Show the image and the classification results 
subplot(1,4,2); 
imshow(uint8(I31)); 
text(10, 20, char(label1),'Color','white');
title('a = 10%');

% Show the image and the classification results 
subplot(1,4,3);
imshow(uint8(I32)); 
text(10, 20, char(label2),'Color','white')
title('a = 20%');


% Show the image and the classification results 
subplot(1,4,4);
imshow(uint8(I33)); 
text(10, 20, char(label3),'Color','white');
title('a = 30%');
