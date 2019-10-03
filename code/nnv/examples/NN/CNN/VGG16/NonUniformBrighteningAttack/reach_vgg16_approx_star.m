
clc;
clear;


fprintf('\n\n=============================LOAD VGG16 ======================\n');

% Load the trained model 
net = vgg16();


fprintf('\n\n=========CONSTRUCT INPUT SET (AN IMAGESTAR SET) =============\n');
% Read the image to classify 
I0 = imread('peppers.png');

% Adjust size of the image 
sz = net.Layers(1).InputSize; 
I = I0(1:sz(1),1:sz(2),1:sz(3));

n = size(I);
N = n(1)*n(2)*n(3);

I4 = reshape(I, [N,1]);
I5 = I4;
del = 20;

for i=1:N
    if I4(i) >= 255 - del
        I5(i) = 255 - I4(i);
    else
        I5(i) = 0;
    end
end

center = reshape(I4, [n(1), n(2), n(3)]); % center image matrix
basis_mat = reshape(I5, [n(1), n(2), n(3)]); % basis image matrix

C = [1;-1];   % 0% <= alpha <= bv percentage of brightening attack
bv = 0.0002;
d = [bv; 0];
pred_lb = 0;
pred_ub = bv;

V(:,:,:,1) = center;
V(:,:,:,2) = basis_mat;

IS = ImageStar(double(V), C, d, pred_lb, pred_ub);


fprintf('\n\n========= PARSE VGG16 FOR REACHABILITY ANALYSIS ============\n');

net = CNN.parse(net, 'VGG16');

fprintf('\n\n======= DO REACHABILITY ANLAYSIS WITH APPROX-STAR METHOD ======\n');

net.reach(IS, 'approx-star');

