onnxfiles = dir("onnx/*.onnx");

for i=1:length(onnxfiles)
    net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name), "InputDataFormats", "BC");
    % nn = matlab2nnv(net);
end

% TODO: just need to implement a slice layer, this should be a priority

