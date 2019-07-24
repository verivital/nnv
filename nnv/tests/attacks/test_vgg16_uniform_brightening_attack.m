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

I3 = reshape(I, [N,1]);

a1 = 0.1; % 10% brightening attack
a2 = 0.5; % 50% brightening attack
a3 = 0.9; % 90% brightening attack

I4 = zeros(N,1, 'uint8');
for i=1:N
    I4(i) = 255 - I3(i);
end


I31 = I3 + a1*I4;
I32 = I3 + a2*I4;
I33 = I3 + a3*I4;

I31 = reshape(I31, [n(1) n(2) n(3)]);
I32 = reshape(I32, [n(1) n(2) n(3)]);
I33 = reshape(I33, [n(1) n(2) n(3)]);

figure;
imshow(I31);
title("Uniformly 20% brightenning attacked image");
figure;
imshow(I32);
title("Uniformly 50% brightenning attacked image");
figure;
imshow(I33);
title("Uniformly 90% brightenning attacked image");


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
