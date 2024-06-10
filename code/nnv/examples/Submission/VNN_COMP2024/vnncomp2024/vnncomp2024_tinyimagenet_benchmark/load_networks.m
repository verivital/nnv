% Load tinyImagenet benchmarks

netPath = "onnx/";

nns = dir(netPath);

for i=1:length(nns)
    if endsWith(nns(i).name, ".onnx")
        disp(nns(i).name);
        onnxNet = importNetworkFromONNX(netPath + nns(i).name, "InputDataFormats","BCSS", "OutputDataFormats","BC");
    end
end

% One network only, it also seems fine