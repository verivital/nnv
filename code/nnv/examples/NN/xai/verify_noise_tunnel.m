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