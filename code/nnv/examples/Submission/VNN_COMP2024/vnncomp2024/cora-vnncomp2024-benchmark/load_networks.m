% Load CORA benchmarks

netPath = "nns/";

nns = dir(netPath);

for i=1:length(nns)
    if endsWith(nns(i).name, ".onnx")
        disp(nns(i).name);
        onnxNet = importNetworkFromONNX(netPath + nns(i).name, "InputDataFormats","BC", "OutputDataFormats","BC");
    end
end

% All good, simple as expected