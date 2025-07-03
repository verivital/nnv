onnxfile = "onnx/model.onnx";

net = importNetworkFromONNX(onnxfile, "InputDataFormats", "BC");
nn = matlab2nnv(net);


