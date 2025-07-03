onnxfiles = dir("onnx/*.onnx");

for i=1:length(onnxfiles)
    if contains(onnxfiles(i).name, "mnist")
        net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name), "InputDataFormats", "BTC");
    elseif contains(onnxfiles(i).name, "oval")
        net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name), "InputDataFormats", "BCSS");
    else
        net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name), "InputDataFormats", "BC");
    end
    nn = matlab2nnv(net);
    % [ ~, name , ext ] = fileparts( onnxfiles(i).name );
    % rmdir("+"+name, 's');
end