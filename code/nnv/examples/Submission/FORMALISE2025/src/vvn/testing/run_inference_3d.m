modelName = "C3D_zoom_out_3d_32x32.onnx";
netonnx = importONNXNetwork("../../models/" + modelName, "InputDataFormats", "TBCSS", "OutputDataFormats", "BC");
net = matlab2nnv(netonnx);
net.OutputSize = 10;
disp("Finished loading model: " + modelName);

% Load data
data = readNPY('../../data/3D/ZoomOut/mnistvideo_32x32_test_data_seq.npy');
labels = readNPY('../../data/3D/ZoomOut/mnistvideo_32x32_test_labels_seq.npy');

tc = 0;

for i=1:length(data)
    s = data(i,:,:,:,:);
    s = permute(s, [1, 3, 2, 4, 5]);
    s = squeeze(s);
    l = labels(i);

    outputs = net.evaluate(s);
    [~, P] = max(outputs);

    if l+1 == P
        tc = tc + 1;
    end
end

fprintf("Total correct: %d \n", tc);
