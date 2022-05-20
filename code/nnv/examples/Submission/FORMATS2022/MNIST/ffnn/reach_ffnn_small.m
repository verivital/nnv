function reach_ffnn_small(pix,numT,noise,XTest,YTest,cora,perturbation)
%% Reachability analysis of an image classification ODE_FFNN (MNIST - small)
% Architecture of first ffnn mnist model:
%  - Inputs = 784 (Flatten images to 1D, original 28x28)
%  - Outputs = 10 (One hot vector)
%  - layer1: relu (64)
%  - layer2: relu (10)
%  - layer3: ODEBlock{
%     -linear (10)
%     }
%  - layer4: linear (10)

%% Section 1. odeblock with NNV reachability

%% Part 1. Loading and constructing the NeuralODE

% Load network parameters
file_path = '../networks/odeffnn_mnist_small.mat';
load(file_path); % Load neuralODe parameters 
% Contruct NeuralODE
layer1 = LayerS(Wb{1},Wb{2}','poslin');
layer2 = LayerS(Wb{3},Wb{4}','poslin');
% ODEBlock only linear layers
% Convert in form of a linear ODE model
states = 10;
outputs = 10;
w1 = Wb{5};
b1 = Wb{6}';
Aout = w1;
Bout = b1;
Cout = eye(states);
D = zeros(outputs,1);
tfinal = 1;
numSteps = 20;
odeblock = LinearODE(Aout,Bout,Cout,D,tfinal,numSteps);
reachStep = tfinal/numSteps;
% Output layers 
layer4 = LayerS(Wb{7},Wb{8}','purelin');
odelayer = ODEblockLayer(odeblock,tfinal,reachStep,false);
neuralLayers = {layer1, layer2, odelayer, layer4};
neuralode = NeuralODE(neuralLayers);

%% Part 2. Load data and prepare experiments

noise = noise*255; % noise 1/10 of max pixel value
% pixels_attack = randi([1 784],1,pix);
pixels_attack = randperm(784,pix);
pred = zeros(numT,1);
time = zeros(numT,1);
for i=1:numT
    img_flat = XTest(:,:,:,i);
    img_flat = reshape(img_flat', [1 784])';
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
    inpS = Star(lb,ub);
    %% Part 3. Reachability 
    t = tic;
    Rode = neuralode.reach(inpS);
    time(i) = toc(t);
    [rv,~] = neuralode.checkRobust(Rode,YTest(i));
    pred(i) = rv;

end

%% Part 3. Robustness results
pred_acc = pred == 1;
timeT = sum(time);
sum_acc = sum(pred_acc);
rob = sum_acc;
acc = sum_acc/numT;
disp(' ');
disp(' ');
disp(' ');
disp('========== Robustness results (NNV) ===========')
disp(' ');
disp('Total images evaluated: '+string(numT));
disp('ATTACK: Random noise, max value = '+string(noise));
disp('Network robust to '+string(sum_acc)+' images.');
disp('Total time to evaluate ' + string(numT) + ' images: ' + string(sum(time)) + ' seconds');
disp('Average time per image: ' + string(sum(time)/numT));

save("ffnn_small_nnv_"+string(perturbation)+"_"+string(noise)+".mat",'rob','pred','timeT','pix','numT','noise');

%% Section 2. ODEblock with CORA reachability
if cora
    %% Part 1. Loading and constructing the NeuralODE

    % Convert in form of a linear ODE model
    odeblockC = LinearODE_cora(Aout,Bout,Cout,D,reachStep,tfinal); % (Non)linear ODE plant
    odelayer = ODEblockLayer(odeblockC,tfinal,reachStep,false);
    neuralLayers = {layer1, layer2, odelayer, layer4};
    neuralode = NeuralODE(neuralLayers);

    %% Part 2. Load data and prepare experiments

    pred = zeros(numT,1);
    timeC = zeros(numT,1);
    for i=1:numT
        img_flat = XTest(:,:,:,i);
        img_flat = reshape(img_flat', [1 784])';
        lb = img_flat;
        ub = img_flat;
        for j=pixels_attack
            ub(j) = min(255, ub(j)+rand*noise);
            lb(j) = max(0, lb(j)-rand*noise);
        end
        % Normalize input (input already normalized)
        lb = lb./255;
        ub = ub./255;
        inpS = Star(lb,ub);
        %% Part 3. Reachability
        t = tic;
        Rode = neuralode.reach(inpS);
        timeC(i) = toc(t);
        [rv,~] = neuralode.checkRobust(Rode,YTest(i));
        pred(i) = rv;
    end


    %% Robustness results

    pred_acc = pred == 1;
    timeT = sum(time);
    sum_acc = sum(pred_acc);
    rob = sum_acc;
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

    save('ffnn_small_cora.mat','rob','pred','timeT','pix','numT','noise');

    % Notify finish
    sound(tan(1:3000));
end
end