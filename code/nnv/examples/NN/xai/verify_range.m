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
% reachOptions.reachMethod = 'exact-star'; % using exact/approx method
reachOptions.reachMethod = 'approx-star';

% Reachability analysis
% R(28,28) = ImageStar;
for i = 1:size(img,1)
    for j = 1:size(img,2)
        lb = img;
        ub = img;
        lb(i,j) = 0;
        ub(i,j) = 1;
        IS = ImageStar(lb,ub);
        t = tic;
        % R{i,j} = net.reach(IS, reachOptions);
        R(i,j) = net.reach(IS, reachOptions);
        toc(t);
    end
end

%% Part 2. Verify XAI

% takes less than a minute to compute this with exact reachability

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
saveas(gcf, "single_pixel_range.png");
