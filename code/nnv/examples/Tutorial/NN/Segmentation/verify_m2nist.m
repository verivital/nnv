% Run a segmentation example using a segnet with transposed convolution

% Load network
net = load("../../../NN/SemanticSegmentation/M2NIST/models/m2nist_75iou_transposedcnn_avgpool.mat");
net = matlab2nnv(net.net);

% Load images
images = load('../../../NN/SemanticSegmentation/M2NIST/m2nist_6484_test_images.mat');
im_data = single(images.im_data);

% Create example input set
Nmax = 50; % maximum allowable number of attacked pixels
de = 0.0001; % disturbance
Nt = 150; % threshold value
% Randomly select 1 images to verify
rng(2);
img_idx = randperm(1000,1);

% Create input set from adversarial perturbation

% Initialize vars
ct = 0; % keep track of pixels modified
flag = 0; % determine when to stop modifying pixels
im = im_data(:,:,img_idx); % select image to evaluate
at_im = im;
% Create brightening attack
for i=1:size(im,1)
    for j=1:size(im,2)
        if im(i,j) > Nt
            at_im(i,j) = 0;
            ct = ct + 1;
            if ct == Nmax
                flag = 1;
                break;
            end
        end
    end
    if flag == 1
        break;
    end
end

% Define input set as ImageStar
dif_im = im - at_im;
noise = -dif_im;
V(:,:,:,1) = im; % center of set
V(:,:,:,2) = noise; % basis vectors
C = [1; -1]; % constraints
d = [1; de-1];
IS = ImageStar(V, C, d, 1-de, 1); % input set
% Label
GrTruth = {im};

% Verify network
reachOptions.reachMethod = 'approx-star';
t = tic;
[riou, rv, rs, n_rb, n_mis, n_unk, n_att, ver_rs, eval_seg_ims] = net.verify_segmentation(IS, GrTruth, reachOptions);
toc(t);

% Visualize results
net.plot_segmentation_output_set(ver_rs{1}, eval_seg_ims{1});



%% Notes
% IoU = Intersection over Union
%  
% IoU = TP / (TP + FP + FN), where
%  
%  - TP = true positive
%  - FP = false positive
%  - FN = false negative
% 
% These metrics are all computed in terms of number of pixels

