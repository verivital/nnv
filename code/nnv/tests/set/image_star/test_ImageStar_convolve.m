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

image1 = ImageStar(IM(:,:,1), LB(:,:,1), UB(:,:,1));

W = [1 -1 0; 0 -1 0; 1 0 1]; % filter
padding = [0 0 0 0];
stride = [1 1];
dilation = [1 1];

conv_image1 = image1.convolve(W, padding, stride, dilation); % convolved image

% for 3 channels image

image = ImageStar(IM, LB, UB); 

W(:,:,1) = [1 -1 0; 0 -1 0; 1 0 1];
W(:,:,2) = [0 1 0; 1 0 1; 0 1 1];
W(:,:,3) = [1 1 0; -1 0 -1; 0 1 1];

conv_image = image.convolve(W, padding, stride, dilation);
