function FilesPath = VCASVerifNeuralNetworks

% Contains the directory of VCAS Verification neural networks data files.
%
% INPUTS
%
% OUTPUTS
%
% FilesPath: cell array containing in each position the directory
% containing the files to be read.

% FilesPath{1} = 'nnv_format/CartPoleSingleReLU.mat';

for i =1:9
FilesPath{i} = [sprintf('nnv_format/nn_vcas_ReLu_%d.mat',i)];
end


end