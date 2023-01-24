function results = verify_tllverify_nnv(onnxF, vnnlibF, reachOpt)
    % Verification of tllverify using NNV
    
    % load network
    loadOpt.InputDataFormat = "BC";
    nn = onnx2nnv(onnxF, loadOpt);
    
    % verify property
    
end

