onnxfiles = dir("onnx/*.onnx");

for i=1:length(onnxfiles)
    
    net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name),     "InputDataFormats","BC","OutputDataFormats","BC");
    nn = matlab2nnv(net);

end