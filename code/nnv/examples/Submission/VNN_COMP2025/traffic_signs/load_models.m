onnxfiles = dir("onnx/*.onnx");

for i=1:length(onnxfiles)
    net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name), "InputDataFormats", "BC");
    % nn = matlab2nnv(net);
end

% could just replace 1 or 2 layers, should not be too complicated
% sign and other simple operations