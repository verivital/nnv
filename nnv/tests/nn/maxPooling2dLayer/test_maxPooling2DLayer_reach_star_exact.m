clear;
clc;
IM(:, :, 1) = [0 0 2 0; 1 2 0 2 ; 0 0 2 2; 0 2 2 2]; % channel 1 input matrix
IM(:, :, 2) = [1 2 2 1; 2 1 2 0; 2 2 2 0; 1 1 1 0]; % channel 2 input matrix
IM(:, :, 3) = [0 0 2 2; 0 2 1 1; 0 2 0 0; 0 2 1 0]; % channel 3 input matrix

LB(:,:,1) = [-0.1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image = ImageStar(IM, LB, UB);

L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

tic;
images = L.reach_star_exact(image);
toc;
