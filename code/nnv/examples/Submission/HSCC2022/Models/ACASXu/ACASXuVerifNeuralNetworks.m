function FilesPath = ACASXuVerifNeuralNetworks

% Contains the directory of ACASXu Verification neural networks data files.
%
% INPUTS
%
% OUTPUTS
%
% FilesPath: cell array containing in each position the directory
% containing the files to be read.

for i =1:5
FilesPath{i} = [sprintf('nnv_format/ACASXU_run2a_%d_1_batch_2000.mat',i)];
end

end