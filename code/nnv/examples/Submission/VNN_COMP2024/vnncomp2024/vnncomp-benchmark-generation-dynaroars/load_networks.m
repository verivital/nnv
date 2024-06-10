% Load dynaroars benchmarks

netPath = "networks/";

nns = dir(netPath);

for i=1:length(nns)
    if endsWith(nns(i).name, ".onnx")
        disp(nns(i).name);
        onnxNet = importNetworkFromONNX(netPath + nns(i).name, "InputDataFormats","BCSS", "OutputDataFormats","BC");
    end
end

% All networks are loaded, no problem here