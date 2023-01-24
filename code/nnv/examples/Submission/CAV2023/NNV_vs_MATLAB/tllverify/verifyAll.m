% function verifyAll()
    % Verify vnnlib properties
    csvFile = "instances.csv";
    opts = detectImportOptions(csvFile);
    opts.Delimiter = ',';
    NNs_props_timeout = readtable(csvFile, opts);
    % only first 24 properties apply to lindex networks, and these
    % alternate (odd -> lindex, even -> lindex_deep)

    % Example to load networks
%     net = importONNXNetwork("onnx/tllBench_n=2_N=M=8_m=1_instance_0_0.onnx", InputDataFormats="BC");
%     loadOpt.InputDataFormat = "BC";
%     nn = onnx2nnv("onnx/tllBench_n=2_N=M=8_m=1_instance_0_0.onnx", loadOpt);
% end

