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


% No need to train or test with perturbed images (I think)

% Create the perturbed images
% training set     -> 25 images per class
% test and val set -> 5 images per class
% These images are added to the original dataset

% Total of 35 images to attack per class (random noise)
nA_train = 25; 
nA_test = 5;
nA_total = 35;
train_attack_idx = randperm(500, nA_total); % only choose from the first 500 images, can reuse indexes for val and test

% Choose randomly which images to attack (same indexes per class)
xattack = zeros(64,64,1,nA_train*6);
yattack = categorical(zeros(nA_train*6,1));
max_pixels = ones(64,64)*255;
min_pixels = zeros(64,64);

% Get extra training data with perturbed images
for i = 1:6 % classes
    iN = (i-1)*2500;
    for j = 1:nA_train % images to attack
        % Select image to apply noise to
        im_idx = train_attack_idx(j)+iN;
        im = Xtrain(:,:,:,im_idx);
        % Apply attack (noise) with a proportion of 0.15 to the image
        noise = (rand(size(Xtrain(:,:,:,im_idx)))-0.5)*25; % equivalent of adding noise of 25 (25 pixel values as images are not normalized)
        idx = rand(size(Xtrain(:,:,:,im_idx))) < 0.15; % less than 15% of the pixels in the image
        im_atk = noise.*idx + im;
        % Ensure upper and lower bounds are within 0 to 255;
        im_atk = min(im_atk, max_pixels);
        im_atk = max(min_pixels, im_atk);
        % Save image with noise
        xattack(:,:,:,(i-1)*6+j) = im_atk;
        yattack((i-1)*6+j) = YData(iN+j);
    end
end

xTest_attack = zeros(64,64,1,nA_test*6);
yTest_attack = categorical(zeros(nA_test*6,1));
xVal_attack = zeros(64,64,1,nA_test*6);
yVal_attack = categorical(zeros(nA_test*6,1));
% Get extra testing and validation data with perturbed images
for i = 1:6 % classes
    iN = (i-1)*500;
    for j = 1:nA_test % images to attack
        % First do test
        % ------------------------
        % Select image to apply noise to
        im_idx = train_attack_idx(j)+iN;
        im = Xtest(:,:,:,im_idx);
        % Apply attack (noise) with a proportion of 0.15 to the image
        noise = (rand(size(Xtest(:,:,:,im_idx)))-0.5)*25; % equivalent of adding noise of 25 (25 pixel values as images are not normalized)
        idx = rand(size(Xtest(:,:,:,im_idx))) < 0.15; % less than 15% of the pixels in the image
        im_atk = noise.*idx + im;
        % Ensure upper and lower bounds are within 0 to 255;
        im_atk = min(im_atk, max_pixels);
        im_atk = max(min_pixels, im_atk);
        % Save image with noise
        xTest_attack(:,:,:,(i-1)*6+j) = im_atk;
        yTest_attack((i-1)*6+j) = Ytest(iN+j);

        % Now do validation data
        % ------------------------
        im_idx = train_attack_idx(j)+iN;
        im = Xval(:,:,:,im_idx);
        % Apply attack (noise) with a proportion of 0.15 to the image
        noise = (rand(size(Xval(:,:,:,im_idx)))-0.5)*25; % equivalent of adding noise of 25 (25 pixel values as images are not normalized)
        idx = rand(size(Xval(:,:,:,im_idx))) < 0.15; % less than 15% of the pixels in the image
        im_atk = noise.*idx + im;
        % Ensure upper and lower bounds are within 0 to 255;
        im_atk = min(im_atk, max_pixels);
        im_atk = max(min_pixels, im_atk);
        % Save image with noise
        xVal_attack(:,:,:,(i-1)*6+j) = im_atk;
        yVal_attack((i-1)*6+j) = Yval(iN+j);
    end
end

% Add attack images to sets
Xtrain = cat(4,Xtrain,xattack);        Ytrain = cat(1, Ytrain, yattack);
Xtest  = cat(4, Xtest, xTest_attack);  Ytest  = cat(1, Ytest, yTest_attack);
Xval   = cat(4, Xval, xVal_attack);    Yval   = cat(1, Yval, yVal_attack);

% Prepare data to pass in training function
data.X = Xtrain;    data.Y = Ytrain;
data.Xtest = Xtest; data.Ytest = Ytest;
data.Xval = Xval;   data.Yval = Yval;


%% Now train all the models

disp("Training models with Jacobian regularization...");

disp("... with glorot initialization...")
medNist_training(data, 'glorot');

disp("... with he initialization...")
medNist_training(data, 'he');

disp("... with narrow-normal initialization...")
medNist_training(data, 'narrow-normal')

disp("Finished training all Jacobian models.")
