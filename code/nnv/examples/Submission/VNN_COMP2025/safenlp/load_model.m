onnxfile1 = "onnx/medical/perturbations_0.onnx";
onnxfile2 = "onnx/ruarobot/perturbations_0.onnx";

net = importNetworkFromONNX(onnxfile1, "InputDataFormats", "BC");
nn = matlab2nnv(net);


net = importNetworkFromONNX(onnxfile2, "InputDataFormats", "BC");
nn = matlab2nnv(net);