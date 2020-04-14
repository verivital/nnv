% this example is from Figure 5 in the website:
% https://ujjwalkarn.me/2016/08/11/intuitive-explanation-convnets/

I = [1 1 1 0 0; 0 1 1 1 0; 0 0 1 1 1; 0 0 1 1 0; 0 1 1 0 0]; % input
W = [1 0 1; 0 1 0; 1 0 1]; % fiter

padding = [0 0 0 0];
stride = [1 1];
dilation = [1 1];


featureMap = Conv2DLayer.compute_featureMap(I, W, padding, stride, dilation);


