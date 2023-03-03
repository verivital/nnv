function FilesPath = CartPoleVerifNeuralNetworks

% Contains the directory of CartPole Verification neural networks data files.
%
% INPUTS
%
% OUTPUTS
%
% FilesPath: cell array containing in each position the directory
% containing the files to be read.

% FilesPath{1} = 'nnv_format/CartPoleSingleReLU.mat';

for i =1:2
FilesPath{i} = [sprintf('nnv_format/nn_cartpole_ReLu_v2_%d.mat',i)];
end


end