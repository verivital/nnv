onnxfiles = dir("onnx/*.onnx");

for i=1:length(onnxfiles)
    net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name), "InputDataFormats", "BCSS");
    nn = matlab2nnv(net);
end