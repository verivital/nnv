function accP = eval_ffnn_large(XTest,YTest)
    %% Test an image classification ODE_FFNN (MNIST)
    % Architecture of first ffnn mnist model:
    %  - Inputs = 784 (Flatten images to 1D, original 28x28)
    %  - Outputs = 10 (One hot vector)
    %  - layer1: relu (64)
    %  - layer2: relu (32)
    %  - layer3: relu (32)
    %  - layer4: relu (32)
    %  - layer5: linear (16)
    %  - layer6: ODEBlock{
    %     -linear (10)
    %     -linear (10)
    %     -linear (10)
    %     -linear (16)
    %     }
    %  - layer7: linear (10)
    
    %% Section 1. odeblock with NNV reachability
    
    %% Part 1. Loading and constructing the NeuralODE
    
    % Load network parameters
    file_path = '../networks/odeffnn_mnist_large.mat';
    load(file_path); % Load neuralODe parameters 
    % Contruct NeuralODE
    layer1 = FullyConnectedLayer(Wb{1}, Wb{2}');
    layer1a = ReluLayer;
    layer2 = FullyConnectedLayer(Wb{3}, Wb{4}');
    layer2a = ReluLayer;
    layer3 = FullyConnectedLayer(Wb{5}, Wb{6}');
    layer3a = ReluLayer;
    layer4 = FullyConnectedLayer(Wb{7}, Wb{8}');
    layer4a = ReluLayer;
    layer5 = FullyConnectedLayer(Wb{9}, Wb{10}');
    % ODEBlock only linear layers
    % Convert in form of a linear ODE model
    states = 16;
    outputs = 16;
    w1 = Wb{11};
    b1 = Wb{12}';
    w2 = Wb{13};
    b2 = Wb{14}';
    w3 = Wb{15};
    b3 = Wb{16}';
    w4 = Wb{17};
    b4 = Wb{18}';
    Aout = w4*w3*w2*w1;
    Bout = w4*(w3*(w2*b1+b2)+b3)+b4;
    Cout = eye(states);
    D = zeros(outputs,1);
    tfinal = 1;
    numSteps = 20;
    odeblock = LinearODE(Aout,Bout,Cout,D,tfinal,numSteps);
    reachStep = tfinal/numSteps;
    odelayer = ODEblockLayer(odeblock,tfinal,reachStep,false);
    % Output layer
    layer7 = FullyConnectedLayer(Wb{19}, Wb{20}');
    % Create neuralode (NN)
    neuralLayers = {layer1, layer1a, layer2, layer2a, layer3, layer3a, layer4, ...
        layer4a, layer5, odelayer, layer7};
    neuralode = NN(neuralLayers);
    
    %% Part 2. Load data and prepare experiments
    numT = length(YTest);
    eval_ode = zeros(numT,1);
    t=tic;
    for i=1:numT
        img_flat = XTest(:,:,:,i);
        img_flat = reshape(img_flat', [1 784])';
        % Normalize input
        img_flat = img_flat./255;
        %% Part 3. Evaluation
        eV = neuralode.evaluate(img_flat);
        [~,eV] = max(eV);
        eval_ode(i) = eV;
    end
    toc(t);
    acc = eval_ode == YTest;
    acc = sum(acc);
    accP = acc/numT;
    disp(' ');
    disp('Correctly classified images: '+string(acc));
    disp('Percentage: '+string(accP));

end

