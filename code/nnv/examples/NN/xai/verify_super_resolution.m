%% Verify the importance of pixels using mnist

%% Part 1. Compute reachability

% Load the model
net = importNetworkFromONNX("super_resolution.onnx", "InputDataFormats","BCSS", "OutputDataFormats","BC");
net = matlab2nnv(net);

% Load data (no download necessary)
digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
    'nndatasets','DigitDataset');
% Images
imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

% Load one image in dataset
[img, fileInfo] = readimage(imds,7010);
target = single(fileInfo.Label); % label = 0 (index 1 for our network)
img = single(img)/255; % change precision

% First, we need to define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using exact/approx method
% does approx method makes sense here?
% use approx for now as exact may compute multiple sets
% at least for now, to show somehing

max_noise = 0.02; % range -> [-0.02, 0.02] (approx 5 pixel color values)

% Reachability analysis
R(28,28) = ImageStar;
for i = 1:size(img,1)
    for j = 1:size(img,2)
        disturbance = zeros(size(img), 'single'); % make a copy 
        disturbance(i,j) = max_noise;
        IS = ImageStar(img, -disturbance, disturbance);
        t = tic;
        R(i,j) = net.reach(IS, reachOptions);
        toc(t);
    end
end

%% Part 2. Verify XAI

% takes less than a minute to compute this with exact reachability

img_scores = net.evaluate(img);

% but what are we verifying exactly? we could compute interval differences
% with original scores (img_scores) to see how things change, but what does
% that mean?

% Let's do the attribution method (baseline) computation, 
% but using the input generation of the noise tunnel

% We cannot use the exact sets we computed, but we can get the ranges for
% all the outputs and use those for the attribution calculation for each
% pixel

ub = zeros(size(img));
lb = zeros(size(img));

for i = 1:size(img,1)
    for j = 1:size(img,2)
        [l,u] = R(i,j).getRange(1,1,target);
        lb(i,j) = l - img_scores(target);
        ub(i,j) = u - img_scores(target);
    end
end

diff = ub - lb;

mapC = hot; % hot map to show the attributions

% Visualize results
figure;
subplot(2,2,1);
imshow(img);
subplot(2,2,2);
imshow(lb, 'Colormap',winter, 'DisplayRange',[min(lb, [], 'all'), max(lb, [], 'all')]);
colorbar;
subplot(2,2,4);
imshow(ub, 'Colormap',winter, 'DisplayRange',[min(ub, [], 'all'), max(ub, [], 'all')]);
colorbar;
subplot(2,2,3);
imshow(diff, 'Colormap',winter, 'DisplayRange',[min(diff, [], 'all'), max(diff, [], 'all')]);
colorbar;

%% Verify 1st method 
% assume something like brightnening or darkening attack for modifying pixel value, e.g. furthest from current value

pixel_val_change = 5/255; % 5 pixel color values for range of pixel

% First, we need to define the reachability options
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star'; % using exact/approx method

% Reachability analysis
R(28,28) = ImageStar;
for i = 1:size(img,1)
    for j = 1:size(img,2)
        % Create set with brightening/darkening pixel
        lb = img; ub = img; % initialize bounds
        if img(i,j) > 122 % darkening
            lb(i,j) = 0;
            ub(i,j) = pixel_val_change;
        else              % brightening
            lb(i,j) = 1-pixel_val_change;
            ub(i,j) = 1;
        end
        % Create ImageStar
        IS = ImageStar(lb,ub);
        % Compute reachable set
        t = tic;
        R(i,j) = net.reach(IS, reachOptions);
        toc(t); % track computation time
    end
end

% Visualize results
img_scores = net.evaluate(img);

ub = zeros(size(img));
lb = zeros(size(img));

for i = 1:size(img,1)
    for j = 1:size(img,2)
        [l,u] = R(i,j).getRange(1,1,target);
        lb(i,j) = l - img_scores(target);
        ub(i,j) = u - img_scores(target);
    end
end

diff = ub - lb;

mapC = hot; % hot map to show the attributions

% Visualize results
figure;
subplot(2,2,1);
imshow(img);
subplot(2,2,2);
imshow(lb, 'Colormap',winter, 'DisplayRange',[min(lb, [], 'all'), max(lb, [], 'all')]);
colorbar;
subplot(2,2,4);
imshow(ub, 'Colormap',winter, 'DisplayRange',[min(ub, [], 'all'), max(ub, [], 'all')]);
colorbar;
subplot(2,2,3);
imshow(diff, 'Colormap',winter, 'DisplayRange',[min(diff, [], 'all'), max(diff, [], 'all')]);
colorbar;


%% Notes
% unsure how to do the noise tunnel one, but we could start with the
% feature ablation

% but the feature ablation just iterates through all pixels, sets them to 0
% and computes some scores... What could we "verify" there?