%
% 2D classification example: MedNIST
% 6 classes: [Abdominal CT, Breast MRI, Chest X-ray, chest CT, Hand, Head CT]
%   labels : [      1          2             3           4       5      6   ]           
%              

%% 1) Load data 

% Download if necessary
dataFolder = "../../../../../../../../data/MedNIST";

% Go through every folder (label) and load all images
categs = dir(dataFolder);
% Initialize vars
XData = zeros(64,64,1,58954); % height = width = 64, 10k imgs per class (greyscale)
YData = zeros(58954, 1);        % except BreatMRI, that has only 8954
% Load images
count = 1;
for i = 3:length(categs)-1
    label = dir(dataFolder + "/"+ string(categs(i).name));
    for k = 3:length(label)
        XData(:, :, :, count) = imread([label(k).folder '/' label(k).name]);
        YData(count) = i-2;
        count = count + 1;
    end
end
YData = categorical(YData);

% Prepare data
rng(0); % fix random seed to select training data
% Let's use 3500 images from each class: 2500 (train), 500 (val), 500 (test)
N = 3500;
idxs = randperm(10000, N); % shuffle indexes for data random
idxs2 = randperm(8954, N); % same but for third class only (less images than others)

% Split data
ntrain = 2500;
ntest = 500; % the other 500 are for validation
% Get all indexes into their corresponding classes
idxs_train = [idxs(1:ntrain), idxs2(1:ntrain)+10000, idxs(1:ntrain)+18954,...
    idxs(1:ntrain)+28954, idxs(1:ntrain)+38954, idxs(1:ntrain)+48954];
idxs_test = [idxs(ntrain+1:ntrain+ntest), idxs2(ntrain+1:ntrain+ntest)+10000, idxs(ntrain+1:ntrain+ntest)+18954,...
    idxs(ntrain+1:ntrain+ntest)+28954, idxs(ntrain+1:ntrain+ntest)+38954, idxs(ntrain+1:ntrain+ntest)+48954];
idxs_val = [idxs(ntrain+ntest+1:end), idxs2(ntrain+ntest+1:end)+10000, idxs(ntrain+ntest+1:end)+18954,...
    idxs(ntrain+ntest+1:end)+28954, idxs(ntrain+ntest+1:end)+38954, idxs(ntrain+ntest+1:end)+48954];

% Allocate data to corresponding variables
Xtrain = XData(:,:,:,idxs_train); Ytrain = YData(idxs_train);
Xtest  = XData(:,:,:, idxs_test); Ytest  = YData(idxs_test);
Xval   = XData(:,:, :, idxs_val); Yval   = YData(idxs_val);

% Prepare data to pass in training function
data.X = Xtrain;    data.Y = Ytrain;
data.Xtest = Xtest; data.Ytest = Ytest;
data.Xval = Xval;   data.Yval = Yval;


%% Now train all the models

disp("Training models with Jacobian regularization...");

% disp("... with glorot initialization...")
% medNist_training(data, 'glorot');
% 
% disp("... with he initialization...")
% medNist_training(data, 'he');

disp("... with narrow-normal initialization...")
medNist_training(data, 'narrow-normal')

disp("Finished training all Jacobian models.")
