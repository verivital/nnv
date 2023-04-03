function reach_cnn_medium(pix,numT,noise,XTest,YTest,perturbation)
    % Reachability analysis of an image classification ODE_CNN (MNIST)
    
    %% 1) Load and construct NN
    
    % Load network parameters
    file_path = 'odecnn_mnist_mid2.mat';
    load(file_path); % Load parameters 

    % layer 1
    w1 = permute(Wb{1},[4 3 2 1]);
    b1 = reshape(Wb{2},[1,1,10]);
    layer1 = Conv2DLayer(w1,b1,[0,0,0,0],[1,1],[1,1]);
    % layer 2
    w2 = reshape(Wb{3},[1 1 10]);
    b2 = reshape(Wb{4},[1 1 10]);
    m2 = reshape(Wb{15},[1 1 10]);
    v2 = reshape(Wb{16},[1 1 10]);
    layer2 = BatchNormalizationLayer('Offset',b2,'Scale',w2,'TrainedMean',m2,'TrainedVariance',v2,'NumChannels',10,'Epsilon',10);
    % layer 3
    layer3 = ReluLayer;
    % layer 4
    w4 = permute(Wb{5},[4 3 2 1]);
    b4 = reshape(Wb{6},[1,1,10]);
    layer4 = Conv2DLayer(w4,b4,[1,1,1,1],[2,2],[1,1]);
    % layer 5
    w5 = reshape(Wb{7},[1 1 10]);
    b5 = reshape(Wb{8},[1 1 10]);
    m5 = reshape(Wb{18},[1 1 10]);
    v5 = reshape(Wb{19},[1 1 10]);
    layer5 = BatchNormalizationLayer('Offset',b5,'Scale',w5,'TrainedMean',m5,'TrainedVariance',v5,'NumChannels',10,'Epsilon',10);
    % layers 6 and 7
    layer6 = ReluLayer;
    layer7 = FlattenLayer;
    layer7.Type = 'nnet.cnn.layer.FlattenLayer';
    % odeblock layer
    states = 1690;
    outputs = 1690;
    w1 = Wb{9};
    b1 = Wb{10}';
    w2 = Wb{11};
    b2 = Wb{12}';
    Aout = w2*w1;
    Bout = w2*b1+b2;
    Cout = eye(states);
    D = zeros(outputs,1);
    tfinal = 1;
    numSteps = 20;
    odeblock = LinearODE(Aout,Bout,Cout,D,tfinal,numSteps);
    reachStep = tfinal/numSteps;
    odelayer = ODEblockLayer(odeblock,1,reachStep);
    % Output layer
    layerout = FullyConnectedLayer(Wb{13}, Wb{14}');
    
    % Create NN
    neuralLayers = {layer1, layer2, layer3, layer4, layer5, layer6, layer7, odelayer, layerout};
    neuralode = NN(neuralLayers);
    
    %% 2) Prepare data and run experiments
    
    % Define attack parameters and initialize variables
    noise = noise*255; % noise 1/10 of max pixel value
    pixels_attack = randperm(784,pix);
    time_ode = zeros(numT,1);
    rob_ode = zeros(numT,1);

    % Create Adversarial attack
    for i=1:numT
        img_flat = XTest(:,:,:,i)';
        lb = img_flat;
        ub = img_flat;
        if strcmp(perturbation,'random')
            for j=pixels_attack
                ub(j) = min(255, ub(j)+noise);
                lb(j) = max(0, lb(j)-noise);
            end
        elseif strcmp(perturbation,'inf')
            for j=pixels_attack
                ub(j) = min(255, ub(j)+noise);
                lb(j) = max(0, lb(j)-noise);
            end
        else
            error('Wrong perturbation type')
        end

        % Normalize input (input already normalized)
        lb = lb./255;
        ub = ub./255;
        inpS = ImageStar(lb,ub);
    
        % Reachability computation
        t = tic;
        Rode = neuralode.reach(inpS);
        time_ode(i) = toc(t);
        rv = neuralode.checkRobust(Rode,YTest(i));
        rob_ode(i) = rv;
    end
    
    % Process results
    rob = sum(rob_ode == 1);
    unk = sum(rob_ode == 0);
    notr = sum(rob_ode == 2);
    timeT = sum(time_ode);
    disp(' ');
    disp('Robust images: '+string(rob));
    disp('Unknown images: '+string(unk));
    disp('Not robust images: '+string(notr));
    disp('Total time = ' + string(timeT));
    % Save results
    save("cnn_medium_nnv_"+string(perturbation)+"_"+string(noise)+".mat",'rob','rob_ode','timeT','pix','numT','noise');

end