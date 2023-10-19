%% First, get the data loaded (same images used for all the models)

% Load data
dataPath = "../../../../../../../data/MNIST/";
filenameImagesTest = 't10k-images-idx3-ubyte.gz';
filenameLabelsTest = 't10k-labels-idx1-ubyte.gz';

% use test data for verification
XData = processImagesMNIST(dataPath + filenameImagesTest);
YData = processLabelsMNIST(dataPath + filenameLabelsTest);
YData = double(YData);

% Get indxs to verify
verInfo = load("acc_results.mat");
XData = XData(:, :, :, verInfo.xVerIdxs);
YData = YData(verInfo.xVerIdxs);

N = length(YData); % total number of images to verify

% 3) Adversarial attack (L_inf)
epsilon = 3/255; 
ub_max = ones(28,28);
lb_min = zeros(28,28);
nR = 100;
xRand = cell(N);

for i=1:N
    img = XData(:,:,:,i);
    lb = img - epsilon;
    lb = max(lb, lb_min); % ensure no negative values
    ub = img + epsilon;
    ub = min(ub, ub_max); % ensure no values > 255 (max pixel value)
    set_idx = (dd-1)*N+i;
    lb = reshape(lb, [28*28, 1]);
    ub = reshape(ub, [28*28, 1]);
    xB = Box(lb, ub); % lb, ub must be vectors
    xImg = xB.sample(nR-2);
    xImg = [lb, ub, xImg];
    xRand{i} = xImg;
end

%% Then verify all data with each model (5*3*3 = 45 models)

% Iterate trhough every folder and subfolder and analyze all the models
path = pwd;
folders = dir(path);
% Skip the first two that appear in every folder and subfolder as those
% correspond to (".", and "..")

% Go into every folder of and analyze each model
for r = 4:6 % iterate through regularizers (3)
    sub_path = [path, filesep, folders(r).name, filesep];
    inits_path = dir(sub_path);
    for i = 3:length(inits_path) % go through all initializations (3 x 3)
        if inits_path(i).isdir
            temp_path = [sub_path, inits_path(i).name, filesep, 'models', filesep];
            models_path = dir([temp_path, '*.mat']);
            for m = 1:length(models_path) % go through all models ( 5 x 3 x 3 )
                netpath = [temp_path, models_path(m).name];
                mnist_falsify_model(netpath, xRand, YData);
            end
        end
    end
end

% 45 different files should be generated under the directory "BenchmarkGen/results_falsify/"

