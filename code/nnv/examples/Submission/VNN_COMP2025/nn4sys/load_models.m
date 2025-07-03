onnxfiles = dir("onnx/*.onnx");

for i=1:length(onnxfiles)
    if contains(onnxfiles(i).name, "lindex")
        net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name), "InputDataFormats", "BC");
    else
        % net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name), "InputDataFormats", "BTC");
        warning("Not supported yet")
    end
    nn = matlab2nnv(net);
end

% Notes: it is actually not that hard, but we'll need to implement some
% slice operations and what now. We'll need to create the networks kinda
% manually, but it is doable. It could be an initial "trial" and then add
% full support. If not by deadline, let's do it right after.

