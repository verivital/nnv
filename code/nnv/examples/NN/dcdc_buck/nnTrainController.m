% Prerequisites: First run buck_parameters.m
% This script generates NN controller after training 
k = 8; % Neurons/layer
S = [k,k,k,k,k,k,k]; % Layers
pid_net = feedforwardnet(S);
for i = 1:length(S)
    pid_net.layers{i}.transferFcn = 'poslin'; % hidden layers poslin, tansig
end

pid_net.layers{length(S)+1}.transferFcn = 'purelin'; % output layer
pid_net.inputs{1}.processFcns = {};
pid_net.outputs{length(S)+1}.processFcns = {};
% config
pid_net = configure(pid_net,input_data,output_data);
pid_net.plotFcns = {'plotperform','plottrainstate',...
    'ploterrhist','plotregression'};
pid_net.trainFcn = 'trainlm';

[x_tot,xi_tot,ai_tot,t_tot] = ...
            preparets(pid_net,input_data,output_data);
pid_net.trainParam.epochs = 800;
pid_net.trainParam.min_grad = 1e-10;
[pid_net,tr] = train(pid_net,x_tot,t_tot,xi_tot,ai_tot);
%generate nn controller 
gensim(pid_net)

%store weights to network
weights(1) = pid_net.IW(1,1);
for i = 1:length(S)
    weights(i+1) = pid_net.LW(i+1,i);
end
bias = pid_net.b;

network.weights = weights;
network.bias = bias;
save('network.mat', 'network');