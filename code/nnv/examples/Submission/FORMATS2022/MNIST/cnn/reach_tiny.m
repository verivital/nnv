function reach_tiny(pix,numT,noise,XTest,YTest,cora,perturbation)
%% Reachability analysis of an image classification ODE_FFNN (MNIST)
% Architecture of first ffnn mnist model:
%  - Inputs = 28x28 images 
%  - Outputs = 10 (One hot vector)
%  - layer1: convolutional (1,4,3,1)
%  - layer2: batchnorm (4)
%  - layer3: relu (no weights)
%  - layer4: convolutional (4,4,3,1)
%  - layer5: batchnorm (4)
%  - layer6: relu (no weights)
%  - layer7: flatten()
%  - layer8: ODEBlock{
%     -linear (10)
%     -linear (676)
%     }
%  - layer9: linear (10)

%% Section 1. odeblock with NNV reachability

%% Part 1. Loading and constructing the NeuralODE

% Load network parameters
file_path = '../networks/odecnn_mnist_tiny.mat';
load(file_path); % Load neuralODe parameters 
% Contruct NeuralODE
% w1 = reshape(Wb{1},[3 3 1 16]);
w1 = permute(Wb{1},[4 3 2 1]);
b1 = reshape(Wb{2},[1,1,4]);
layer1 = Conv2DLayer(w1,b1,[0,0,0,0],[1,1],[1,1]);
w2 = reshape(Wb{3},[1 1 4]);
b2 = reshape(Wb{4},[1 1 4]);
m2 = reshape(Wb{15},[1 1 4]);
v2 = reshape(Wb{16},[1 1 4]);
layer2 = BatchNormalizationLayer('Offset',b2,'Scale',w2,'TrainedMean',m2,'TrainedVariance',v2,'NumChannels',4,'Epsilon',4);
layer3 = ReluLayer;
% w4 = reshape(Wb{5},[4 4 16 16]);
w4 = permute(Wb{5},[4 3 2 1]);
b4 = reshape(Wb{6},[1,1,4]);
layer4 = Conv2DLayer(w4,b4,[1,1,1,1],[2,2],[1,1]);
w5 = reshape(Wb{7},[1 1 4]);
b5 = reshape(Wb{8},[1 1 4]);
m5 = reshape(Wb{18},[1 1 4]);
v5 = reshape(Wb{19},[1 1 4]);
layer5 = BatchNormalizationLayer('Offset',b5,'Scale',w5,'TrainedMean',m5,'TrainedVariance',v5,'NumChannels',4,'Epsilon',4);
layer6 = ReluLayer;
layer7 = FlattenLayer;
layer7.Type = 'nnet.cnn.layer.FlattenLayer';
% layer7.Type = 'nnet.keras.layer.FlattenCStyleLayer';
layers = {layer1, layer2, layer3, layer4, layer5, layer6, layer7};
net1 = CNN('cnn',layers,1,1); % neural network controller
% ODEBlock only linear layers
% Convert in form of a linear ODE model
states = 676;
outputs = 676;
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
% Output layers 
layerout = LayerS(Wb{13},Wb{14}','purelin');
layer_out = FFNNS(layerout);

odelayer = ODEblockLayer(odeblock,1,reachStep,false);
neuralLayers = {layer1, layer2, layer3, layer4, layer5, layer6, layer7, odelayer, layerout};
neuralode = NeuralODE(neuralLayers);

%% Part 2. Load data and prepare experiments

noise = noise*255; % noise 1/10 of max pixel value
% pixels_attack = randi([28 28],1,pix);
pixels_attack = randperm(784,pix);
pred = zeros(numT,1);
time = zeros(numT,1);
pred_ode = zeros(numT,1);
time_ode = zeros(numT,1);
% InpSS = []; % Array of input sets
rob_ode = zeros(numT,1);
for i=1:numT
%     img_flat = double(XTest(:,:,:,i));
%     img_flat = extractdata(img_flat)';
%     img_flat = reshape(img_flat', [1 784])';
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
%     img_inp = img_flat./255;
    %% Part 3. Reachability and Simulation
    t = tic;
    Rode = neuralode.reach(inpS);
    time_ode(i) = toc(t);
    [rv,~] = neuralode.checkRobust(Rode,YTest(i));
    rob_ode(i) = rv;
end

rob = sum(rob_ode == 1);
unk = sum(rob_ode == 0);
notr = sum(rob_ode == 2);
timeT = sum(time_ode);
disp(' ');
disp('Robust images: '+string(rob));
disp('Unknown images: '+string(unk));
disp('Not robust images: '+string(notr));
disp('Total time = ' + string(timeT));

save("cnn_tiny_nnv_"+string(perturbation)+"_"+string(noise)+".mat",'rob','rob_ode','timeT','pix','numT','noise');

% Notify finish
sound(tan(1:3000));


%% Section 2. ODEblock with CORA reachability
if cora
    %% Part 1. Loading and constructing the NeuralODE

    % Convert in form of a linear ODE model
    % C = eye(states); % Want to get both of the outputs from NeuralODE
    odeblockC = LinearODE_cora(Aout,Bout,Cout,D,reachStep,tfinal); % (Non)linear ODE plant 

    %% Part 2. Load data and prepare experiments

    pred = zeros(numT,1);
    timeC = zeros(numT,1);
    for i=1:numT
%         img_flat = double(XTest(:,:,:,i));
%         img_flat = extractdata(img_flat)';
    %     img_flat = reshape(img_flat', [1 784])';
        img_flat = XTest(:,:,:,i);
        lb = img_flat;
        ub = img_flat;
        for j=pixels_attack
            ub(j) = min(255, ub(j)+rand*noise);
            lb(j) = max(0, lb(j)-rand*noise);
        end
        % Normalize input (input already normalized)
        lb = lb./255;
        ub = ub./255;
        inpS = ImageStar(lb,ub);
        img_inp = img_flat./255;
        %% Part 3. Reachability and Simulation
        % Divide reachability into steps
        U = Star(0,0);
        t = tic;
        R1 = net1.reach(inpS,'approx-star');
        R1 = R1.toStar;
        R2 = odeblockC.stepReachStar(R1,U);
        R3 = layer_out.reach(R2,'approx-star');
        time(i) = toc(t);
        [lb_out,ub_out] = R3.getRanges;
        [maxO,max_idx] = max(ub_out);
        pred(i) = max_idx;
    end


    %% Robustness results

    pred_acc = pred == YTest(1:numT);
    sum_acc = sum(pred_acc);
    acc = sum_acc/numT;
    disp(' ');
    disp(' ');
    disp('========== Robustness results (CORA) ===========')
    disp(' ');
    disp('Total images evaluated: '+string(numT));
    disp('ATTACK: Random noise, max value = '+string(noise));
    disp('Network robust to '+string(sum_acc)+' images.');
    disp('Total time to evaluate ' + string(numT) + ' images: ' + string(sum(time)) + ' seconds');
    disp('Average time per image: ' + string(sum(time)/numT));
end

% Notify finish
sound(tan(1:3000));
end