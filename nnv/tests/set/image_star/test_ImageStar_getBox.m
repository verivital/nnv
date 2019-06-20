% test constructor for ImageStar class
IM(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1]; % center image channel 1
IM(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1]; % center image channel 2
IM(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0]; % center image channel 3

LB(:,:,1) = [-0.1 -0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; % attack on pixel (1,1) and (1,2)
LB(:,:,2) = [-0.1 -0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]; 
LB(:,:,3) = LB(:,:,2);

UB(:,:,1) = [0.1 0.2 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,2) = [0.1 0.15 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB(:,:,3) = UB(:,:,2);

image = ImageStar(IM, LB, UB);
tic;
[image_lb, image_ub] = image.getBox;
toc;
display(image_lb);
display(image_ub);

% NOTE***: using parallel computation may not help
% NOTE***: turning on the parallel pool adds noise, the computation time
% increase even we don't use parallel pool for computation