function medNist_training(data, initializerP)
    %
    % 2D classification example: MedNIST
    % 6 classes: [Abdominal CT, Breast MRI, Chest X-ray, chest CT, Hand, Head CT]
    %   labels : [      1          2             3           4       5      6   ]           
    % 
    % Training code based on:
    % https://www.mathworks.com/help/deeplearning/ug/train-deep-learning-network-with-jacobian-regularization.html
    % 
    
    % 1) Get data
    % divide data again
    Xtrain = data.X;     Ytrain = data.Y;
    Xtest = dlarray(data.Xtest, "SSCB");  Ytest = double(data.Ytest);
    Xval  = dlarray(data.Xval, "SSCB");   Yval = double(data.Yval);

    %% 2) Create model
    
    % Strat with a small model and see how it performs
    
    % Define layers
    if ~strcmp(initializerP, "narrow-normal")
        layers = [imageInputLayer([64 64 1], Mean=mean(Xtrain,4))
            convolution2dLayer(3, 3,'Stride',1,"WeightsInitializer",initializerP, 'BiasInitializer', 'zeros')
            batchNormalizationLayer
            reluLayer
            averagePooling2dLayer(2,'Stride',2,'Padding',[0 0 0 1])
            flattenLayer
            fullyConnectedLayer(6, "WeightsInitializer",initializerP, 'BiasInitializer', 'zeros')
            softmaxLayer];
    else
        layers = [imageInputLayer([64 64 1], Mean=mean(Xtrain,4))
            convolution2dLayer(3, 3,'Stride',1,"WeightsInitializer",initializerP, 'BiasInitializer',initializerP)
            batchNormalizationLayer
            reluLayer
            averagePooling2dLayer(2,'Stride',2,'Padding',[0 0 0 1])
            flattenLayer
            fullyConnectedLayer(6, "WeightsInitializer",initializerP, 'BiasInitializer',initializerP)
            softmaxLayer];
    end
    
    net = dlnetwork(layers);
    
    
    %% Specify training options
    
    jacobianRegularizationCoefficient = 0.01; % MATLAB example = 1
    learningRate = 0.01;
    momentum = 0.9;
    
    numEpochs = 5;
    miniBatchSize = 128;
    xN = length(Ytrain); % number of samples used in training

    
    %% 3) Train model
    
    seeds = [0,1,2,3,4];
    
    for sd = seeds
        
        rng(sd);
        acc_val = 0; % initialize validation accuracy
        iteration = 0;
        velocity = [];
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
                T = Ytrain(batch_idxs(batch_range));
                % One-hot encode labels
                T = onehotencode(T',1);
             
                % Evaluate the model loss, gradients and the network state using
                % dlfeval and the modelLoss function listed at the end of the example.
                [~, gradTotalLoss, state] = dlfeval(@modelLoss, net, X, T, ...
                    miniBatchSize, jacobianRegularizationCoefficient);
                net.State = state;
                
                % Update the network parameters.
                [net, velocity] = sgdmupdate(net,gradTotalLoss,velocity,learningRate,momentum);
    
                y_val = predict(net,Xval);
                [~, y_val] = max(extractdata(y_val), [], 1);
                acc_val_iter = sum(y_val' == Yval)/numel(Yval);
                if acc_val < acc_val_iter
                    best_net = net; % save best trained network
                    acc_val = acc_val_iter;
                end
                
                % Plot the training progress.
                if ~mod(iteration,50)
                    D = duration(0,0,toc(start),Format="hh:mm:ss");
                    disp("Training with Jacobian regularization" + ", Epoch: " + epoch + ", Elapsed: " + string(D) + ", Accuracy: " + string(acc_val_iter));
                end
                
            end
        end
        
        
        %% 4) Test model
        
        net = best_net;
    
        y = predict(net,Xtest);
        [~, y] = max(extractdata(y), [], 1);
        accuracy = sum(y' == Ytest)/numel(Ytest);
        disp("Test accuracy = " + string(accuracy));
        
        % Save model
        save(initializerP + "/models/model_jacobian_"+initializerP+"_"+string(sd)+".mat", 'net', 'accuracy');
    
    end
end

%% Helper functions
function [totalLoss, gradTotalLoss, state] = modelLoss(net, X, T, miniBatchSize, jacobianRegularizationCoefficient)
    
    % Find prediction and loss.
    [Z,state] = forward(net, X);
    loss = crossentropy(Z, T);
    
    numClasses = size(Z,1);
    numProjections = 1;
    regularizationTerm = 0;
    
    % Compute Jacobian term and its gradient.
    for i = 1:numProjections
       
        % Sample a matrix whose elements are drawn from the standard Normal distribution.
        rndarray = randn(numClasses, miniBatchSize);
        
        % Normalize the columns of the random matrix.
        rndarray = normc(rndarray);
        
        % Compute the vector-vector product.
        vectorproduct = rndarray(:)' * Z(:);  
       
        % Compute the gradient of the vector-vector product. Since another
        % derivative will be taken, set EnableHigherDerivatives to true.
        vectorJacobianTerm = dlgradient(vectorproduct, X, EnableHigherDerivatives=true);
        
        % Multiply by necessary constants to obtain approximation of the
        % Frobenius norm of the Jacobian.
        regularizationTerm = regularizationTerm + numClasses*sum(vectorJacobianTerm.^2,"all") /(numProjections*miniBatchSize);
    end
    totalLoss = loss + jacobianRegularizationCoefficient/2 * regularizationTerm;
    
    % Calculate the gradient of the loss.
    gradTotalLoss = dlgradient(totalLoss, net.Learnables);
end

