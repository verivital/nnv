function accP = eval_tiny(XTest,YTest)
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

%% Part 1. Loading and constructing the NeuralODE

% Load network parameters
file_path = '../../../../Python/neuralODE_examples/mnist/odecnn_mnist_tiny.mat';
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

odelayer = ODEblockLayer(odeblock,1,reachStep,false);
neuralLayers = {layer1, layer2, layer3, layer4, layer5, layer6, layer7, odelayer, layerout};
neuralode = NeuralODE(neuralLayers);

%% Part 2. Load data and prepare experiments
numT = length(YTest);
eval_ode = zeros(numT,1);
t = tic;
for i=1:numT
    img_t = XTest(:,:,:,i)';
    % Normalize input
    img_t = img_t./255;
    %% Part 3. Evaluation
    eV = neuralode.evaluate(img_t);
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

