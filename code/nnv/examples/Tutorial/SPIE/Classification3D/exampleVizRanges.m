%% Visualization example using the organCmnist

% This data set has 11 possible classification labels
% 1)  heart
% 2)  left lung
% 3)  right lung
% 4)  liver
% 5)  spleen
% 6)  pancreas
% 7)  left kidney
% 8)  right kidney
% 9)  bladder
% 10) left femoral head
% 11) right femoral head

% references
% https://ieee-dataport.org/documents/annotations-body-organ-localization-based-miccai-lits-dataset#files
% https://medmnist.com/
% From the medmnist journal publication:
%   "Organ{A,C,S}MNIST. Te Organ{A,C,S}MNIST is based on 3D computed tomography (CT) images from
%    Liver Tumor Segmentation Benchmark (LiTS)30. Tey are renamed from OrganMNIST_{Axial,Coronal,Sagittal}
%    (in MedMNIST v19) for simplicity. We use bounding-box annotations of 11 body organs from another study31 to
%    obtain the organ labels. Hounsfeld-Unit (HU) of the 3D images are transformed into gray-scale with an abdominal window. 
%    We crop 2D images from the center slices of the 3D bounding boxes in axial/coronal/sagittal views
%    (planes). Te only diferences of Organ{A,C,S}MNIST are the views. Te images are resized into 1×28×28 to
%    perform multi-class classifcation of 11 body organs. 115 and 16 CT scans from the source training set are used
%    as training and validation set, respectively. Te 70 CT scans from the source test set are treated as the test set.


%% Robustness verification of a NN (L infinity adversarial attack)
%  if f(x) = y, then forall x' in X s.t. ||x - x'||_{\infty} <= eps,
%  then f(x') = y = f(x)

% Load network 
nn_model = load('models/model_organmnist3d.mat');

% Create NNV model
net = matlab2nnv(nn_model.net);

dataset = "data/organmnist3d.mat";

% Load data
load(dataset); % only the test set

% data to verify (test set)
test_images = permute(test_images, [2 3 4 1]);
test_labels = test_labels + 1;

% adversarial attack
epsilon = 1; % {epsilon} color values

% select image N to verify
N = 5; % ?
img = single(test_images(:,:,:,N));
target = test_labels(N);
numClasses = 11;

% volshow(img); % to visualize input

% Create input set
% ensure the values are within the valid range for pixels ([0 255])
lb_min = zeros(size(img)); % minimum allowed values for lower bound is 0
ub_max = 255*ones(size(img)); % maximum allowed values for upper bound is 255
lb_clip = max((img-epsilon),lb_min);
ub_clip = min((img+epsilon), ub_max);
IS = VolumeStar(lb_clip, ub_clip); % this is the input set we will use

% Let's evaluate the image and the lower and upper bounds to ensure these
% are correctly classified

% Evaluate input image
Y_outputs = net.evaluate(img); 
[~, yPred] = max(Y_outputs); % (expected: y = 1)

% Evaluate lower and upper bounds
LB_outputs = net.evaluate(lb_clip);
[~, LB_Pred] = max(LB_outputs); % (expected: y = 1)
UB_outputs = net.evaluate(ub_clip);
[~, UB_Pred] = max(UB_outputs); % (expected: y = 1)

% Now, we can do the verification process of this image w/ L_inf attack

% The easiest way to do it is using the verify_robustness function

% First, we need to define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'relax-star-area'; % using relax method
reachOptions.relaxFactor = 1;

% Verification
t = tic;
res_approx = net.verify_robustness(IS, reachOptions, target);

if res_approx == 1
    disp("Neural network is verified to be robust!")
else
    disp("Unknown result")
end

toc(t);


%% Let's visualize the ranges for every possible output

R = net.reachSet{end};

[lb_out, ub_out] = R.getRanges;
lb_out = squeeze(lb_out);
ub_out = squeeze(ub_out);

mid_range = (lb_out + ub_out)/2;

range_size = ub_out - mid_range;

x = [1 2 3 4 5 6 7 8 9 10 11];

figure;
errorbar(x, mid_range, range_size, '.');
hold on;
xlim([0.5 11.5]);
scatter(x,Y_outputs, 'x', 'MarkerEdgeColor', 'r');


%% Notes
% The ranges obtained are an overappoximation (projection) of the 
% true ranges of the computed star sets

