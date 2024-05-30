% Load LSNC benchmarks

netPath = "onnx/";

nns = dir(netPath);

for i=1:length(nns)
    if endsWith(nns(i).name, ".onnx")
        disp(nns(i).name);
        onnxNet = importNetworkFromONNX(netPath + nns(i).name, "InputDataFormats","BC", "OutputDataFormats","BC");
    end
end

% We get these warnings... Not great, probably customized part of the
% network function during training
% Warning: The ONNX file uses IR version 8, while the highest fully-supported IR is version 7. 
% Warning: The ONNX file uses Opset version 16, while the fully-supported version is 6~14. The imported network may differ from the ONNX
% network. 