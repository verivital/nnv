% Load data
% classes = airplane, automobile, bird, cat, deer, dog, frog, horse, ship, and truck.
load("cifar_data_example.mat");
numClasses = 10;

% load network
load("net_cifar3.mat");

epsilon = 1/255;

lb = single(x) - epsilon;
ub = single(x) + epsilon;

IS = ImageStar(lb,ub);
target = 7;

reachOptions.reachMethod = "cp-star"; % not actually needed
reachOptions.train_device = "cpu"; % default = pu, comment out to run on gpu

result = verify_robustness_cp(net,IS,reachOptions,target, numClasses);