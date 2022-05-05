function accP = eval_ffnn(XTest,YTest)
%% Test an image classification ODE_FFNN (MNIST)
% Architecture of first ffnn mnist model:
%  - Inputs = 784 (Flatten images to 1D, original 28x28)
%  - Outputs = 10 (One hot vector)
%  - layer1: relu (64)
%  - layer2: relu (32)
%  - layer3: linear (16)
%  - layer4: ODEBlock{
%     -linear (10)
%     -linear (16)
%     }
%  - layer5: linear (10)

%% Section 1. odeblock with NNV reachability

%% Part 1. Loading and constructing the NeuralODE

% Load network parameters
file_path = '../../../../Python/neuralODE_examples/mnist/odeffnn_mnist.mat';
load(file_path); % Load neuralODe parameters 
% Contruct NeuralODE
layer1 = LayerS(Wb{1},Wb{2}','poslin');
layer2 = LayerS(Wb{3},Wb{4}','poslin');
layer3 = LayerS(Wb{5},Wb{6}','purelin');
% ODEBlock only linear layers
% Convert in form of a linear ODE model
states = 16;
outputs = 16;
w1 = Wb{7};
b1 = Wb{8}';
w2 = Wb{9};
b2 = Wb{10}';
Aout = w2*w1;
Bout = w2*b1+b2;
Cout = eye(states);
D = zeros(outputs,1);
tfinal = 1;
numSteps = 20;
odeblock = LinearODE(Aout,Bout,Cout,D,tfinal,numSteps);
reachStep = tfinal/numSteps;
% Output layers 
layer5 = LayerS(Wb{11},Wb{12}','purelin');
odelayer = ODEblockLayer(odeblock,tfinal,reachStep,false);
neuralLayers = {layer1, layer2, layer3, odelayer, layer5};
neuralode = NeuralODE(neuralLayers);

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

