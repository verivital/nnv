function NNStruct = CreateNNStruct(files, format, reachmethod)

% Loads the NN and creates a struct containing NNs in the correct format.
%
% INPUTS
%
% files: cell array containing the path where NN info is kept.
% format: class that transforms the content in files into the correct
%          format. 
% reachmethod: class containing reachability method.
%
% OUTPUT
%
% NNStruct: struct containing the NNs in the format given by adapter.

numofNN = length(files);

for j = 1:numofNN
    load(files{j});
    n = length(b);
    Layers = [];
    try
        if act_fcns
            for i=1:n
                L = reachmethod(W{i}, b{i},act_fcns{i});
                Layers = [Layers L];
            end
        end
    catch
        for i=1:n
            if i == n
                L = reachmethod(W{i}, b{i},'purelin');
                Layers = [Layers L];
            else
                L = reachmethod(W{i}, b{i},'poslin');
                Layers = [Layers L];
            end
        end
    end
    controller = format(Layers);
    name = sprintf('n%d',j);
    NNStruct.(name) = controller;
  
end
end