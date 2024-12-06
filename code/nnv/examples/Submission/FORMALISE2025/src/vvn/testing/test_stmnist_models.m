%% Load things
% Number of frames
nF = 4;

% Load data
data = readNPY(sprintf("../../data/STMNIST/test/stmnistvideo_%df_test_data_seq.npy", nF));
labels = readNPY(sprintf("../../data/STMNIST/test/stmnistvideo_%df_test_labels_seq.npy", nF));

% Preprocessing
% from [B D C H W] to [B D H W C]
reshaped_data = permute(data, [1, 3, 4, 5, 2]);
datacopy = reshaped_data(:,:,:,:,:);

% Experimental variables
numClasses = 10;

% Load the model
modelName = sprintf("stmnist_%df.onnx", nF);
netonnx = importONNXNetwork("../../models/" + modelName, "InputDataFormats", "BCSSS", "OutputDataFormats", "BC");
net = matlab2nnv(netonnx);
net.OutputSize = numClasses;
disp("Finished loading model: " + modelName);


%% Make predictions on test set to check that we are verifying correct
outputLabels = zeros(length(datacopy), 1);
correct = 0;
total = length(datacopy);

for i=1:total
    s = datacopy(i,:,:,:,:);
    s = squeeze(s);
    l = labels(i) + 1;

    output = net.evaluate(s);
    [~, P] = max(output);

    if P == l
        correct = correct + 1;
    end
end

acc = correct / total;
fprintf("All done! \n")