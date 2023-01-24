% function verify_lindex()
    
    %% Load networks
    % NNV (to show we can directly load them from onnx)
    % lindex
    nn1 = onnx2nnv("onnx/lindex.onnx");
    % lindex_deep
    nn2 = onnx2nnv("onnx/lindex_deep.onnx");
    % MATLAB
    net1 = importONNXNetwork("onnx/lindex.onnx", "InputDataFormats","BC");
    Layers = net1.Layers(1:end-1);
    net1 = dlnetwork(Layers);
    net2 = importONNXNetwork("onnx/lindex_deep.onnx", "InputDataFormats","BC");
    Layers = net2.Layers(1:end-1);
    net2 = dlnetwork(Layers);

    % Verify vnnlib properties
    csvFile = "instances.csv";
    opts = detectImportOptions(csvFile);
    opts.Delimiter = ',';
    NNs_props_timeout = readtable(csvFile, opts);
    % only first 24 properties apply to lindex networks, and these
    % alternate (odd -> lindex, even -> lindex_deep)


    
% end

