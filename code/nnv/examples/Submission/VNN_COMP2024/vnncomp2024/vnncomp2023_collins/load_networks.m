% Load Collins benchmarks

netPath = "onnx/";

nns = dir(netPath);

for i=1:length(nns)
    if endsWith(nns(i).name, ".onnx")
        disp(nns(i).name);
        onnxNet = importNetworkFromONNX(netPath + nns(i).name, "InputDataFormats","BCSS", "OutputDataFormats","BCT");
    end
end

% 1) OutputDataFormats (empty, or BCT)
% Probably not, we can skip this one. 3 reshape to reshape layers, not
% standard, lots of parameters in these custom layers...
% Could be accomplished with some manual conversion for these...
% 2) OutputDataFormats (TBC)
% Only generated 1 custom layer (reshape_to_concat), this may be doable
% But still looks like a very challenging benchmark