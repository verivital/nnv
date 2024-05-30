% Load LinearizeNN benchmarks

netPath = "onnx/AllInOne.onnx";

onnxNet = importNetworkFromONNX(netPath, "InputDataFormats","BC", "OutputDataFormats","BC");

% All good as well, only one network here