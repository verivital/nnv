%% Image input set examples
% Show examples from an image input set

% Load data (no download necessary)
digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
    'nndatasets','DigitDataset');
% Images
imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

% Load one image in dataset
[img, fileInfo] = readimage(imds,10);
target = double(fileInfo.Label); % label = 0 (index 1 for our network)
img = double(img); % convert to double

% Create input set
disturbance = 5; % 5 pixel color values
lb_min = zeros(size(img)); % minimum allowed values for lower bound is 0
ub_max = 255*ones(size(img)); % maximum allowed values for upper bound is 255
lb_clip = max((img-disturbance),lb_min);
ub_clip = min((img+disturbance), ub_max);
IS = ImageStar(lb_clip, ub_clip); % this is the input set we will use

% Sample some images from the input set
rand_images = IS.sample(3); % get 3 random images from the input set

% Visualize image;
figure;
subplot(4,2,[1 2]);
imshow(img, [0 255]);
title("Original");
subplot(4,2,3);
imshow(rand_images{1}, [0 255]);
title("Sample 1");
subplot(4,2,5);
imshow(rand_images{2}, [0 255]);
title("Sample 2");
subplot(4,2,7);
imshow(rand_images{3}, [0 255]);
title("Sample 3");
subplot(4,2,4);
imshow(img - rand_images{1}, [-5 5]);
title("Diff 1");
colorbar;
subplot(4,2,6);
imshow(img - rand_images{2}, [-5 5]);
title("Diff 2");
colorbar;
subplot(4,2,8);
imshow(img - rand_images{3}, [-5 5]);
title("Diff 3");
colorbar;



