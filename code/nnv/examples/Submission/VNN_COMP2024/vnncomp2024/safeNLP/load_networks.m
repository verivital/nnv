% Load safeNLP benchmarks

netPath = "onnx/";

nns = dir(netPath);

for i=1:length(nns)
    if endsWith(nns(i).name, ".onnx")
        disp(nns(i).name);
        onnxNet = importNetworkFromONNX(netPath + nns(i).name, "InputDataFormats","BC", "OutputDataFormats","BC");
    end
end

% Networks seem very simple... Are we missing something?