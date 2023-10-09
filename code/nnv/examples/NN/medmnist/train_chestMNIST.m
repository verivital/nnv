%% ChestMNIST is a multi-label (14) binary classification problem

t = tic; % track total time for the training

dataset = "data/mat_files/chestmnist.mat";

% Load data
load(dataset);

% data information
numClasses = size(val_labels,2);
imgSize = size(val_images); % size of the images 
imgSize = imgSize(2:end); % first dimension corresponds to number of images

% pre-process data
% convert data to single precision (desired precision for matlab training)
% convert labels to categorical
% permute dimensions of images

% train
train_images = single(train_images);
train_images = permute(train_images, [2 3 4 1]);
Xtrain = dlarray(train_images, "SSCB");
Ytrain = single(train_labels);

% test
test_images = single(test_images);
test_images = permute(test_images, [2 3 4 1]);
Xtest = dlarray(test_images, "SSCB");
Ytest = single(test_labels);

% validation
val_images = single(val_images);
val_images = permute(val_images, [2 3 4 1]);
Xval = dlarray(val_images, "SSCB");
Yval = single(val_labels)';

% Create the neural network model
layers = [
    imageInputLayer(imgSize) % image size = [28 28 1] (Height x Width x Channels)
    
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(numClasses) 
    sigmoidLayer];

net = dlnetwork(layers);

% Implement training
rng(0); % set random seed for reproducibility
numEpochs = 10;
net = train_multilabel(net, Xtrain, Ytrain, Xval, Yval, numEpochs);

% Test network (accuracy)
YPred = predict(net,test_images);

accuracy = sum(YPred == test_labels, 'all')/numel(test_labels);
disp("Test accuracy = "+string(accuracy));

% Save model
name = split(dataset, filesep);
name = name{end};
name = split(name, '.');
name = name{1};
save("models/model_"+string(name)+".mat", 'net', 'accuracy');

toc(t);



%% Helper functions

% Train the model
function best_net = train_multilabel(net, Xtrain, Ytrain, Xval, Yval, numEpochs)

    miniBatchSize = 128;
    iteration = 0;
    acc_val = 0;
    xN = length(Ytrain); % number of samples used in training

    % Initialize parameter for adam optimizer
    trailingAvg = [];
    trailingAvgSq = [];
    learnRate = 0.01;
    gradientDecayFactor = 0.5;
    squaredGradientDecayFactor = 0.999;
    
    start = tic;
    % Loop over epochs.
    for epoch = 1:numEpochs
        
        % Reset and shuffle mini-batch queue.
        batch_idxs = randperm(xN, xN);
        iter_per_batch = floor(length(Ytrain)/miniBatchSize);
        
        for iter = 1:iter_per_batch
            iteration = iteration + 1;
            
            % Read mini-batch of data.
            batch_range = ((iter-1) * miniBatchSize + 1) : ((iter-1) * miniBatchSize + miniBatchSize) ;
            X = Xtrain(:,:,:, batch_idxs(batch_range));
            X = dlarray(X, "SSCB");
            T = Ytrain(batch_idxs(batch_range),:)';

            % If training on a GPU, then convert data to gpuArray.
            if canUseGPU
                X = gpuArray(X);
            end
         
            % Evaluate the model loss, gradients and the network state using
            % dlfeval and the modelLoss function listed at the end of the example.
            [gradients, loss] = dlfeval(@modelLoss, net, X, T);
            
            % Update the network parameters using the Adam optimizer.
            [net.Learnables, trailingAvg, trailingAvgSq] = adamupdate(net.Learnables,gradients, ...
                trailingAvg,trailingAvgSq,iteration,learnRate, ...
                gradientDecayFactor,squaredGradientDecayFactor);

            y_val = extractdata(predict(net,Xval));
            acc_val_iter = sum(y_val == Yval, 'all')/numel(Yval);
            if acc_val < acc_val_iter
                best_net = net; % save best trained network
                acc_val = acc_val_iter;
            end
            
            % Plot the training progress.
            if ~mod(iteration, 100)
                D = duration(0,0,toc(start),Format="hh:mm:ss");
                disp("Training progress... " + " Epoch: " + epoch + ", Elapsed: " + string(D) + ...
                    ", Accuracy: " + string(acc_val_iter) + ", Loss: " + string(loss));
            end
            
        end
    end
    
end


% Compute the loss for every batch
function [gradients, loss] = modelLoss(net, X, T)
    
    % Find prediction and loss.
    Z = forward(net, X);

    % compute loss
    loss = crossentropy(Z, T, ClassificationMode="multilabel");
    
    % Calculate the gradient of the loss.
    gradients = dlgradient(loss, net.Learnables);
end
