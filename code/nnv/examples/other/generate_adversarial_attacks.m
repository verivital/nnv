

%prerequistes and shared variables
I = imread('peppers.png');
I1 = im2double(I);
figure;
imshow(I1);
title("Original Image");

n = size(I1); % size of image
N = numel(I1); % number of elements in image

%used throughout
I_RESHAPED = reshape(I1, [N 1]); %SERENE_CODE

%used by brightening attacks
a1 = 0.2; % 10% brightening attack
a2 = 0.5; % 50% brightening attack
a3 = 0.9; % 90% brightening attack

%used by plot tests
I4 = I_RESHAPED;
I5 = I_RESHAPED;
delta = 0.5;

for i=1:N
    if I4(i) >= 1 - delta
        I5(i) = 1 - I4(i);
    else
        I5(i) = 0;
    end
end

%used by FGSM
a = 0.0;
b = 1;
r = (b-a).*rand(N,1) + a;


%% example 1: brightening attack AI2
% brightenning attack in AI2 (nonuniform brightening attack)

I2 = I_RESHAPED;
delta = 0.5;
for i=1:N
    if I2(i) >= 1-delta
        I2(i) = 1.0;
    end
end

I2 = reshape(I2, [n(1) n(2) n(3)]);
figure;
imshow(I2);
title("Nonuniform brightenning attacked image");


%% example 2: uniform brightening attack
% our new brightening attack on all pixel(uniform brightening attack)
% xi = xi^0 + a*(1-xi^0); a is percentage of attack from 0 - 100 %

I3 = I_RESHAPED;
I31 = I3 + a1*(ones(N,1) - I3);
I31 = reshape(I31, [n(1) n(2) n(3)]);

figure;
imshow(I31);
title("Uniformly 20% brightenning attacked image");


%% example 3: uniform brightening attack
% our new brightening attack on all pixel(uniform brightening attack)
% xi = xi^0 + a*(1-xi^0); a is percentage of attack from 0 - 100 %

I3 = I_RESHAPED;
I32 = I3 + a2*(ones(N,1) - I3);
I32 = reshape(I32, [n(1) n(2) n(3)]);

figure;
imshow(I32);
title("Uniformly 50% brightenning attacked image");


%% example 4: uniform brightening attack
% our new brightening attack on all pixel(uniform brightening attack)
% xi = xi^0 + a*(1-xi^0); a is percentage of attack from 0 - 100 %

I3 = I_RESHAPED;
I33 = I3 + a3*(ones(N,1) - I3);
I33 = reshape(I33, [n(1) n(2) n(3)]);

figure;
imshow(I33);
title("Uniformly 90% brightenning attacked image");

%% example 5: imagestar
% uniform brightening attack

V(:,:,:,1) = I1; % center image matrix
V(:,:,:,2) = ones(n(1),n(2),n(3))-I1; % basis image matrix
C = [1;-1];
d = [1; 0];
pred_lb = 0;
pred_ub = 1;

IS1 = ImageStar(V, C, d, pred_lb, pred_ub);


%% example 6: imagestar
% nonuniform (partially) brightening attack

V(:,:,:,1) = reshape(I4, [n(1), n(2), n(3)]); % center image matrix
V(:,:,:,2) = reshape(I5, [n(1), n(2), n(3)]); % basis image matrix

C = [1;-1];
d = [1; 0];
pred_lb = 0;
pred_ub = 1;

IS2 = ImageStar(V, C, d, pred_lb, pred_ub);

%% example 7: nonuniform brightening attack
% plot nonuniform brightening attack

I61 = I4 + a1*I5;
I61 = reshape(I61, [n(1), n(2), n(3)]);

figure;
imshow(I61);
title("Nonuniformly 20% brightenning attacked image");

%% example 8: nonuniform brightening attack
% plot nonuniform brightening attack

I62 = I4 + a2*I5;
I62 = reshape(I62, [n(1), n(2), n(3)]);

figure;
imshow(I62);
title("Nonuniformly 50% brightenning attacked image");

%% example 9: nonuniform brightening attack
% plot nonuniform brightening attack

I63 = I4 + a3*I5; 
I63 = reshape(I63, [n(1), n(2), n(3)]);

figure;
imshow(I63);
title("Nonuniformly 90% brightenning attacked image");


%% example 10: FGSM attack

delta = 0.1;
I31 = I_RESHAPED + delta*r;
I31 = reshape(I31, [n(1), n(2), n(3)]);

figure; 
imshow(I31);
title("FGSM attack with delta = 0.1");


%% example 11: FGSM attack

delta = 0.2;
I32 = I_RESHAPED + delta*r;
I32 = reshape(I32, [n(1), n(2), n(3)]);

figure; 
imshow(I32);
title("FGSM attack with delta = 0.2");


%% example 12: FGSM attack

delta = 0.3;
I33 = I_RESHAPED + delta*r;
I33 = reshape(I33, [n(1), n(2), n(3)]);

figure; 
imshow(I33);
title("FGSM attack with delta = 0.3");

