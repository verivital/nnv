%% 1) Load data 

%% 1) Load data

dataPath = "../../../../../../../../data/MNIST/";
filenameImagesTrain = 'train-images-idx3-ubyte.gz';
filenameLabelsTrain = 'train-labels-idx1-ubyte.gz';

XData = processImagesMNIST(dataPath + filenameImagesTrain);
YData = processLabelsMNIST(dataPath + filenameLabelsTrain);

% Prepare data
rng(0); % fix random seed to select training data
% Let's use 3500 images from each class: 2500 (train), 500 (val), 500 (test)
N = 3500;
idxs = randperm(6000, N);
ntrain = 2500;
ntest = 500;
% Get all indexes into their corresponding classes
idxs_train = [];
idxs_test = [];
idxs_val = [];
for i=1:10 % 10 classes
    class_train = idxs(1:ntrain) + 6000*(i-1);
    idxs_train = [idxs_train, class_train];
    class_test = idxs(ntrain+1:ntrain+ntest) + 6000*(i-1);
    idxs_test = [idxs_test, class_test];
    class_val = idxs(ntrain+ntest+1:end) + 6000*(i-1);
    idxs_val= [idxs_val, class_val];
end

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

disp("... with glorot initialization...")
mnist_training(data, 'glorot');

disp("... with he initialization...")
mnist_training(data, 'he');

disp("... with narrow-normal initialization...")
mnist_training(data, 'narrow-normal')

disp("Finished training all Jacobian models.")
