%% Function for training of a 3D medmnist 
% Code based on MathWorks example
% https://www.mathworks.com/help/vision/ug/train-classification-network-to-classify-object-in-3-d-point-cloud.html

function train_medmnist3d(dataset)

    disp("Training model on dataset: "+string(dataset));
    
    t = tic; % track total time for the training
    
    % Load data
    load(dataset);

    % data information
    numClasses = length(unique(val_labels)); % number of classes in dataset
    imgSize = size(val_images); % size of the images 
    imgSize = imgSize(2:end); % first dimension corresponds to number of images
    
    % pre-process data
    % convert data to single precision (desired precision for matlab training)
    % convert labels to categorical
    % permute dimensions of images

    % train
    train_images = single(train_images);
    train_images = permute(train_images, [2 3 4 5 1]);
    train_labels = categorical(train_labels + 1);

    % test
    test_images = single(test_images);
    test_images = permute(test_images, [2 3 4 5 1]);
    test_labels = categorical(test_labels + 1);

    % validation
    val_images = single(val_images);
    val_images = permute(val_images, [2 3 4 5 1]);
    val_labels = categorical(val_labels + 1);
    
    
    % Create the neural network model
    layers = [
        image3dInputLayer(imgSize) % image size = [28 28 28] (Height x Width x Channels)
        
        convolution3dLayer(5,28,'Padding','same')
        batchNormalizationLayer
        reluLayer
        averagePooling3dLayer(2,'Stride',2)

        convolution3dLayer(3,28,'Padding','same')
        batchNormalizationLayer
        reluLayer
        averagePooling3dLayer(2,'Stride',2)

        flattenLayer

        fullyConnectedLayer(128)
        reluLayer
        dropoutLayer(0.5)
        
        fullyConnectedLayer(numClasses) % 10 = number of classes
        softmaxLayer
        classificationLayer];
    
    % Training options
    options = trainingOptions('adam', ...
        'InitialLearnRate',0.001, ...
        'MaxEpochs',30, ...
        'MiniBatchSize', 32,...
        'Shuffle','every-epoch', ...
        'ValidationData',{val_images, val_labels}, ...
        'ValidationFrequency',25, ...
        'ValidationPatience', 5,...
        'OutputNetwork','best-validation-loss',...
        'Verbose',true);
    
    % Train network
    [net, info] = trainNetwork(train_images, train_labels, layers, options);
    
    % Test network (accuracy)
    YPred = classify(net,test_images);
    
    accuracy = sum(YPred == test_labels)/numel(test_labels);
    disp("Test accuracy = "+string(accuracy));
    
    % Save model
    name = split(dataset, filesep);
    name = name{end};
    name = split(name, '.');
    name = name{1};
    save("models/model_"+string(name)+".mat", 'net', 'accuracy', 'info');
    
    toc(t);

end
