clear
clc
close all
%addpath(fullfile(matlabroot,'examples','mpc','main'));
global mdl Ts T G_ego t_gap D_default v_set amin_ego amax_ego
global x0_lead v0_lead x0_ego v0_ego seed
mdl = 'mpcACCsystem';
open_system(mdl)

% Generate training data
[inputs_cell, output_cell] = generateOnetrace();
inputs_test = inputs_cell;
output_test = output_cell;
for i = 1:20
    [input1, output1] = generateOnetrace();
    inputs_cell = catsamples(inputs_cell, input1, 'pad');
    output_cell = catsamples(output_cell, output1, 'pad');
end


%% Construct NN controller

S = [10,10,10];
mrac_net = feedforwardnet(S);
for i = 1:length(S)
    mrac_net.layers{i}.transferFcn = 'poslin';
end

mrac_net.layers{length(S)+1}.transferFcn = 'purelin';
mrac_net.inputs{1}.processFcns = {};
mrac_net.outputs{length(S)+1}.processFcns = {};

% mrac_net.layerConnect = [0 0 0 1;1 0 0 0;0 1 0 0;0 0 1 0];
%  mrac_net.layerWeights{1,4}.delays = 1:2;
 
mrac_net = configure(mrac_net,inputs_cell,output_cell);
mrac_net.plotFcns = {'plotperform','plottrainstate',...
    'ploterrhist','plotregression','plotresponse'};
mrac_net.trainFcn = 'trainlm';

[x_tot,xi_tot,ai_tot,t_tot] = ...
            preparets(mrac_net,inputs_cell,output_cell);
mrac_net.trainParam.epochs = 500;
mrac_net.trainParam.min_grad = 1e-10;
[mrac_net,tr] = train(mrac_net,x_tot,t_tot,xi_tot,ai_tot);

testoutseq = mrac_net(inputs_test);
testout = cell2mat(testoutseq);
figure
plot(testout,'r');
hold on 
plot(cell2mat(output_test),'b')
hold off

%generate simlink nn
%gensim(mrac_net)

% store weights to network
% weights(1) = mrac_net.IW(1,1);
% for i = 1:length(S)
%     weights(i+1) = mrac_net.LW(i+1,i);
% end
% bias = mrac_net.b;
% 
% network.weights = weights;
% network.bias = bias;
% save network

%% Remove example file folder from MATLAB path, and close Simulink model.
rmpath(fullfile(matlabroot,'examples','mpc','main'));