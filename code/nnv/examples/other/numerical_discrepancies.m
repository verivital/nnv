%% There are some discrepancies that we can observe here

% Load network 
mnist_model = load('../Tutorial/NN/MNIST/mnist_model_fc.mat');

% Create NNV model
net1 = matlab2nnv(mnist_model.net);
net2 = matlab2nnv(mnist_model.net);
net3 = matlab2nnv(mnist_model.net);
net4 = matlab2nnv(mnist_model.net);

% Load data (no download necessary)
digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
    'nndatasets','DigitDataset');
% Images
imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

% Load first image in dataset
[img, fileInfo] = readimage(imds,600);
img = single(img); % change precision

% Create input set
ones_ = ones(size(img), 'single');
disturbance = 2 .* ones_; % one pixel value (images are not normalized, they get normalized in the ImageInputLayer)

% However, we need to ensure the values are within the valid range for pixels ([0 255])
lb_min = zeros(size(img)); % minimum allowed values for lower bound is 0
ub_max = 255*ones(size(img)); % maximum allowed values for upper bound is 255
lb_clip = max((img-disturbance),lb_min);
ub_clip = min((img+disturbance), ub_max);
IS = ImageStar(lb_clip, ub_clip); % this is the input set we will use

% First, we need to define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approxiate method


%% Let's compute reach sets to evaluate the differences

R1 = net1.reach(IS, reachOptions);
y1 = net1.evaluate(img);

reachOptions.device = 'gpu';
IS = ImageStar(lb_clip, ub_clip); % this is the input set we will use
R2 = net2.reach(IS, reachOptions);
y2 = net2.evaluate(gpuArray(img));


%% Now on double precision
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approximate method

IS = ImageStar(lb_clip, ub_clip); % this is the input set we will use
IS = IS.changeVarsPrecision('double');
net3 = net3.changeParamsPrecision('double');
R3 = net3.reach(IS, reachOptions);
y3 = net3.evaluate(double(img));

reachOptions.device = 'gpu';
ImageStar(lb_clip, ub_clip); % this is the input set we will use
IS = IS.changeVarsPrecision('double');
net4.changeParamsPrecision('double');
R4 = net4.reach(IS, reachOptions);
y4 = net4.evaluate(gpuArray(double(img)));

%% What do the output ranges look like for all of them?

[lb1, ub1] = R1.estimateRanges;
[lb2, ub2] = R2.estimateRanges;
[lb3, ub3] = R3.estimateRanges;
[lb4, ub4] = R4.estimateRanges;

outputRanges = [lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4];
outputRanges = squeeze(outputRanges)'

%% where did it go wrong?
%% Layer 1
r1 = net1.reachSet{1};
r2 = net2.reachSet{1};
r3 = net3.reachSet{1};
r4 = net4.reachSet{1};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer1 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% % layer1 = squeeze(layer1)';
% 
%% Layer 2
r1 = net1.reachSet{2};
r2 = net2.reachSet{2};
r3 = net3.reachSet{2};
r4 = net4.reachSet{2};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer2 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% 
%% Layer 3
r1 = net1.reachSet{3};
r2 = net2.reachSet{3};
r3 = net3.reachSet{3};
r4 = net4.reachSet{3};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer3 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% 
%% Layer 4
r1 = net1.reachSet{4};
r2 = net2.reachSet{4};
r3 = net3.reachSet{4};
r4 = net4.reachSet{4};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer4 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% 
%% Layer 5
r1 = net1.reachSet{5};
r2 = net2.reachSet{5};
r3 = net3.reachSet{5};
r4 = net4.reachSet{5};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer5 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% 
% 
%% Layer 6
r1 = net1.reachSet{6};
r2 = net2.reachSet{6};
r3 = net3.reachSet{6};
r4 = net4.reachSet{6};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer6 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% 
% 
%% Layer 7
r1 = net1.reachSet{7};
r2 = net2.reachSet{7};
r3 = net3.reachSet{7};
r4 = net4.reachSet{7};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer7 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% 
% 
%% Layer 8
r1 = net1.reachSet{8};
r2 = net2.reachSet{8};
r3 = net3.reachSet{8};
r4 = net4.reachSet{8};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer8 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% 
% 
%% Layer 9
r1 = net1.reachSet{9};
r2 = net2.reachSet{9};
r3 = net3.reachSet{9};
r4 = net4.reachSet{9};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer9 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% 
%% Layer 10
r1 = net1.reachSet{10};
r2 = net2.reachSet{10};
r3 = net3.reachSet{10};
r4 = net4.reachSet{10};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer10 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
% 
%% Layer 11
r1 = net1.reachSet{11};
r2 = net2.reachSet{11};
r3 = net3.reachSet{11};
r4 = net4.reachSet{11};
% 
% [lb1, ub1] = r1.estimateRanges;
% [lb2, ub2] = r2.estimateRanges;
% [lb3, ub3] = r3.estimateRanges;
% [lb4, ub4] = r4.estimateRanges;
% 
% layer11 = {lb1; lb2; lb3; lb4; ub1; ub2; ub3; ub4};
