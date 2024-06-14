%% Robustness verification of a NN (L infinity adversarial attack)
%  if f(x) = y, then forall x' in X s.t. ||x - x'||_{\infty} <= eps,
%  then f(x') = y = f(x)

% Load network 
model = load('gtsrb_model.mat');

% Create NNV model
net = matlab2nnv(model.net);

% Load data
gtsrb_path = [nnvroot(), filesep, 'data', filesep, 'GTSRB', filesep];
imds = imageDatastore(gtsrb_path, 'IncludeSubfolders',true,'LabelSource','foldernames');

inputSize = [30 29];
imds.ReadFcn = @(loc)imresize(imread(loc),inputSize);

% Load one image in dataset
[img, fileInfo] = readimage(imds,10);

% Visualize image 
figure;
imshow(img);

% Get image info
target = single(fileInfo.Label); % target label
img = single(img); % change precision


% Create input set

% One way to define it is using original image +- disturbance (L_inf epsilon)
ones_ = ones(size(img), 'single');
disturbance = 1 .* ones_; % one pixel value (images are not normalized, they get normalized in the ImageInputLayer)
I = ImageStar(img, -disturbance, disturbance);

% Can also define it with just lower and upper bounds
I2 = ImageStar(img-disturbance, img+disturbance);

% However, we need to ensure the values are within the valid range for pixels ([0 255])
lb_min = zeros(size(img)); % minimum allowed values for lower bound is 0
ub_max = 255*ones(size(img)); % maximum allowed values for upper bound is 255
lb_clip = max((img-disturbance),lb_min);
ub_clip = min((img+disturbance), ub_max);
IS = ImageStar(lb_clip, ub_clip); % this is the input set we will use


% Let's evaluate the image and the lower and upper bounds to ensure these
% are correctly classified

% Evaluate input image
Y_outputs = net.evaluate(img); 
[~, yPred] = max(Y_outputs); % (expected: y = target)

% Evaluate lower and upper bounds
LB_outputs = net.evaluate(lb_clip);
[~, LB_Pred] = max(LB_outputs); % (expected: y = target)
UB_outputs = net.evaluate(ub_clip);
[~, UB_Pred] = max(UB_outputs); % (expected: y = target)

% Now, we can do the verification process of this image w/ L_inf attack

% The easiest way to do it is using the verify_robustness function

% First, we need to define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using approxiate method

% Verification
t = tic;
res_approx = net.verify_robustness(IS, reachOptions, target);

if res_approx == 1
    disp("Neural network is verified to be robust!")
elseif any([yPred, LB_Pred, UB_Pred] ~= target)
    disp("Neural network is not robust!");
else
    disp("Unknown result")
end

toc(t);


%% Let's visualize the ranges for every possible output

% Get output reachable set
R = net.reachSet{end};

% Get (overapproximate) ranges for each output index
[lb_out, ub_out] = R.getRanges;
lb_out = squeeze(lb_out);
ub_out = squeeze(ub_out);

% Get middle point for each output and range sizes
mid_range = (lb_out + ub_out)/2;
range_size = ub_out - mid_range;

% Label for x-axis
x = 1:net.OutputSize;

% Visualize set ranges and evaluation points
figure;
errorbar(x, mid_range, range_size, '.');
hold on;
scatter(x,Y_outputs, 'x', 'MarkerEdgeColor', 'r');
