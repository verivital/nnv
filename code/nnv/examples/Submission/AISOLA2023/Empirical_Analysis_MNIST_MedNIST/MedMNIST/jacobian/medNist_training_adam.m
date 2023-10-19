function medNist_training_adam(data, initializerP)
    %
    % 2D classification example: MedNIST
    % 6 classes: [Abdominal CT, Breast MRI, Chest X-ray, chest CT, Hand, Head CT]
    %   labels : [      1          2             3           4       5      6   ]           
    % 
    % Training code based on:
    % https://www.mathworks.com/help/deeplearning/ug/train-deep-learning-network-with-jacobian-regularization.html
    % 
    %
    %
%
% I think there is a bug in the adamupdate code (line 92 in this script):
% Unrecognized property 'Value' for class 'dlarray'.
% 
% Error in deep.internal.recording.containerfeval>iProcessNetwork_Nout_Nin (line 354)
%     protoTable.Value = zeros(height(protoTable),0);
% 
% Error in deep.internal.recording.containerfeval>iDispatch_Nout_Nin (line 194)
%     outputs = iProcessNetwork_Nout_Nin(fun, paramFun, numOut, ...
% 
% Error in deep.internal.recording.containerfeval (line 38)
%     outputs = iDispatch_Nout_Nin(allowNetInput, fun, paramFun, numOut, ...
% 
% Error in deep.internal.networkContainerFixedArgsFun (line 29)
% varargout = deep.internal.recording.containerfeval(...
% 
% Error in adamupdate (line 152)
% [p, avg_g, avg_gsq] = deep.internal.networkContainerFixedArgsFun(func, ...
%
%
% They try to assign some value to a dlarray, which not possible because
% dlarrays do not have those properties
    
    % 1) Get data
    % divide data again
    Xtrain = data.X;     Ytrain = data.Y;
    Xtest = data.Xtest;  Ytest = data.Ytest;
    Xval = data.Xval;    Yval = data.Yval;

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
        layers = [imageInputLayer([64 64 1], Mean=mean(XTrain,4))
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
    
    jacobianRegularizationCoefficient = 0.01; % MATLAB examples = 1
    
    numEpochs = 5;
    miniBatchSize = 32;
    xN = length(Ytrain); % number of samples used in training

    
    %% 3) Train model
    
    seeds = [0,1,2,3,4];
    
    for sd = seeds
        
        rng(sd);
        acc_val = 0; % initialize validation accuracy
        iteration = 0;
        gradientsAvg = [];
        sqGradientsAvg = [];
        start = tic;
        
        % Loop over epochs.
        for epoch = 1:numEpochs
            
            % Reset and shuffle mini-batch queue.
            batch_idxs = randperm(xN, xN);
            iter_per_batch = floor(length(Ytrain)/miniBatchSize);
            
            for iter = 1:iter_per_batch
                iteration = iteration + 1;
                
                % Read mini-batch of data.
                batch_range = ((iter-1) * iter_per_batch + iter) : ((iter-1) * iter_per_batch + miniBatchSize) ;
                X = Xtrain(:,:,:, batch_idxs(batch_range));
                X = dlarray(X, "SSCB");
                T = Ytrain(batch_idxs(batch_range));
                % One-hot encode labels
                T = onehotencode(T',1);
             
                % Evaluate the model loss, gradients and the network state using
                % dlfeval and the modelLoss function listed at the end of the example.
                gradients = dlfeval(@modelLoss, net, X, T, ...
                    miniBatchSize, jacobianRegularizationCoefficient);
                
                % Update the network parameters.
                [net, gradientsAvg, sqGradientsAvg] = adamupdate(net, gradients, gradientsAvg, sqGradientsAvg, iteration);
    
                y_val = classify(net,Xval);
                acc_val_iter = sum(y_val == Yval)/numel(Yval);
                if acc_val < acc_val_iter
                    best_net = net; % save best trained network
                    acc_val = acc_val_iter;
                end
                
                % Plot the training progress.
                D = duration(0,0,toc(start),Format="hh:mm:ss");
                disp("Training with Jacobian regularization" + ", Epoch: " + epoch + ", Elapsed: " + string(D))
                
            end
        end
        
        
        %% 4) Test model
        
        net = best_net;
    
        y = classify(net,Xtest);
        accuracy = sum(y == Ytest)/numel(Ytest);
        disp("Test accuracy = " + string(accuracy));
        
        % Save model
        save("models/model_jacobian_"+initializerP+"_"+string(sd)+".mat", 'net', 'accuracy');
    
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

