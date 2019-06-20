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
[image_lb, image_ub] = image.getBox;
display(image_lb(:,:,1)); % display lower bound of channel 1
display(image_ub(:,:,1)); % display upper bound of channel 1

center = [1;1];
others = [1;3];
channel_ind = 1;
tic; 
image1 = image.isMax(center, others, channel_ind); % check if x[1,1] >= x[1,3]
toc;

tic;
image2 = image.isMax(others, center, channel_ind); % check if x[1,3] >= x[1,1]
toc;

x2 = [2; 2];
x3 = [1;3];
others = [x2 x3];

tic;
image3 = image.isMax(center, others, channel_ind); % check if x[1,1] >= x[2,2] && x[1,1] >= x[1,3]
toc; 