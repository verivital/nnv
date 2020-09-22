%to run this as a test, use results_attacks=runtests('test_attacks')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%prerequistes and shared variables

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

I_RESHAPED = reshape(I, [N,1]);
I_RESHAPED_DOUBLE = double(I_RESHAPED);

beta = 255;
r_int = randi(beta,N,1, 'uint8');
r_double = rand(N, 1);


I4 = zeros(N,1, 'uint8');
for i=1:N
    I4(i) = 255 - I_RESHAPED(i);
end


%____________________________________________________________________
%BELOW: tests originally from test_vgg16_construct_imageStar_inputSet.m


%% test 1: vgg16 FGSM imageStar inputSet
% FGSM attack



del1 = 0.2;
I31 = I_RESHAPED + del1*r_int;

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



% TODO: make subtests for stuff like the clasification instead of just
% checking if it runs?


%____________________________________________________________________
%BELOW: tests originally from test_vgg16_FGSM_attack.m

%% test 2: vgg16 FGSM attack display
% Show the image and the classification results 
figure; 
subplot(1,4,1);
imshow(I);
title('Original image');


%% test 3: vgg16 FGSM attack 10 percent
% Show the image and the classification results 
del1 = 0.1;
I31 = I_RESHAPED_DOUBLE + del1*r_double;
I31 = reshape(I31, [n(1), n(2), n(3)]);
label1 = classify(net, I31);
subplot(1,4,2); 
imshow(uint8(I31)); 
text(10, 20, char(label1),'Color','white');
title('a = 10%');


%% test 4: vgg16 FGSM attack 20 percent
% Show the image and the classification results 
del2 = 0.1593543;
I32 = I_RESHAPED_DOUBLE + del2*r_double;
I32 = reshape(I32, [n(1), n(2), n(3)]);
label2 = classify(net, I32);
subplot(1,4,3);
imshow(uint8(I32)); 
text(10, 20, char(label2),'Color','white')
title('a = 20%');


%% test 5: vgg16 FGSM attack 30 percent
% Show the image and the classification results
del3 = 0.5;
I33 = I_RESHAPED_DOUBLE + del3*r_double;
I33 = reshape(I33, [n(1), n(2), n(3)]);
label3 = classify(net, I33);
subplot(1,4,4);
imshow(uint8(I33)); 
text(10, 20, char(label3),'Color','white');
title('a = 30%');


%____________________________________________________________________
%BELOW: tests originally from test_vgg16_get_output_eachLayer.m


%% test 6: vgg16 classification
% Classify the image using VGG-16 
label = classify(net, I);

% Show the image and the classification results 
figure; 
imshow(I) 
text(10, 20, char(label),'Color','white')



%% test 7: vgg16 FGSM output layer
% FGSM attack
del1 = 0.1;
I31 = I_RESHAPED + del1*r_int;
I31 = reshape(I31, [n(1), n(2), n(3)]);

figure; 
imshow(I31);
title("FGSM attack with del = 0.1");

label1 = classify(net, I31);

% Show the image and the classification results 
figure; 
imshow(I31) 
text(10, 20, char(label1),'Color','white')



%% test 8: vgg16 FGSM output layer
% FGSM attack
del2 = 0.2;
I32 = I_RESHAPED + del2*r_int;
I32 = reshape(I32, [n(1), n(2), n(3)]);

figure; 
imshow(I32);
title("FGSM attack with del = 0.2");

label2 = classify(net, I32);

% Show the image and the classification results 
figure; 
imshow(I32) 
text(10, 20, char(label2),'Color','white')



%% test 9: vgg16 FGSM output layer
% FGSM attack
del3 = 0.3;
I33 = I_RESHAPED + del3*r_int;
I33 = reshape(I33, [n(1), n(2), n(3)]);

figure; 
imshow(I33);
title("FGSM attack with del = 0.3");

label3 = classify(net, I33);

% Show the image and the classification results 
figure; 
imshow(I33) 
text(10, 20, char(label3),'Color','white')

%TODO: figure out if this is just doing a run test, or if we are
%trying to confirm that it's resilient or something
%NOT SURE HOW THIS IS DIFFERENT FROM THE EARLIER TESTS


%____________________________________________________________________
%BELOW: tests originally from test_vgg16_uniform_brightening_attack.m


%% test 10: vgg16 uniform brightening
a1 = 0.1; % 10% brightening attack
I31 = I_RESHAPED + a1*I4;
I31 = reshape(I31, [n(1) n(2) n(3)]);

figure;
imshow(I31);
title("Uniformly 20% brightenning attacked image");
label1 = classify(net, I31);

% Show the image and the classification results 
figure; 
imshow(I31) 
text(10, 20, char(label1),'Color','white')


%% test 11: vgg16 uniform brightening
a2 = 0.5; % 50% brightening attack
I32 = I_RESHAPED + a2*I4;
I32 = reshape(I32, [n(1) n(2) n(3)]);

figure;
imshow(I32);
title("Uniformly 50% brightenning attacked image");
label2 = classify(net, I32);

% Show the image and the classification results 
figure; 
imshow(I32) 
text(10, 20, char(label2),'Color','white')



%% test 12: vgg16 uniform brightening
a3 = 0.9; % 90% brightening attack
I33 = I_RESHAPED + a3*I4;
I33 = reshape(I33, [n(1) n(2) n(3)]);

figure;
imshow(I33);
title("Uniformly 90% brightenning attacked image");
label3 = classify(net, I33);

% Show the image and the classification results 
figure; 
imshow(I33) 
text(10, 20, char(label3),'Color','white')

%TODO: figure out if this is just doing a run test, or if we are
%trying to confirm that it's resilient or something



%____________________________________________________________________
%BELOW: tests originally from test_vgg16_zero_center_normalization.m


%% test 13: vgg16 zero center normalization
del1 = 0.2;
I31 = I_RESHAPED + del1*r_int;
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

%well.....it's done but like.....what happened
