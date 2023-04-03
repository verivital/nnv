function [networks, names2idxs] = load_acasxu_NNs()
%% 1) Load networks (collins benchmarks)
% vnnFolder = "/home/manzand/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";
vnnFolder = "/home/dieman95/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";
benchmarkFolder = "acasxu/onnx/";
listNN = dir(vnnFolder+benchmarkFolder);
loadOptions = struct;
loadOptions.InputDataFormat ='BCSS';
networks = {}; % create a cell array of neural networks
names = [];
idxs = [];
count = 1;
t = tic;
for h = 1:length(listNN) % generlize NN loading options for all benchmarks
    if endsWith(listNN(h).name, ".onnx")
        try
            networks{count} = onnx2nnv(vnnFolder+benchmarkFolder+string(listNN(h).name), loadOptions);
        catch 
            networks{count} = "Error";
        end
        names{count} = listNN(h).name;
        idxs{count} = count;
        count = count + 1;
    end
end
t = toc(t);
names2idxs = containers.Map(names,idxs);
disp("All networks are loaded in " + string(t) + " seconds");
end