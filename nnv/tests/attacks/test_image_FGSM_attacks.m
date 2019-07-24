% shearing attack use shear mapping 
% the information of shear mapping is here: https://en.wikipedia.org/wiki/Shear_mapping

I = imread('peppers.png');
I1 = im2double(I);
figure;
imshow(I1);
title("Original Image");

n = size(I1);
N = n(1)*n(2)*n(3);


% nonuniform (partially) brightening attack
I2 = reshape(I1, [N,1]);

a = 0.0;
b = 1;
r = (b-a).*rand(N,1) + a;
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

