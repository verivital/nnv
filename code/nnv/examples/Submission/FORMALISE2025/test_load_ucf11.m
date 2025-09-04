modelName = sprintf("ucf11_c3d_%df.onnx", 16);
% netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "TBCSS", "OutputDataFormats", "BC");
netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "BCSSS", "OutputDataFormats", "BC");
net = matlab2nnv(netonnx);
net.OutputSize = 11;
disp("Finished loading model: " + modelName);

data = readNPY(sprintf("data/UCF11/ucf11_data_normalized_%df.npy", 16)); % 100 16 112 112 3
labels = readNPY(sprintf("data/UCF11/ucf11_labels_%df.npy", 16));
s = data(1,:,:,:,:);
s = squeeze(s);

output = net.evaluate(s);

% TESTED AND OUTPUTS MATCH
%% 

disp(size(s)); % 8 120 160 3
VS = L_inf_attack(s, 1/255, 8);


%% Helper Functions
function VS = L_inf_attack(x, epsilon, numFrames)
    lb = squeeze(x);
    ub = squeeze(x);

    % Perturb the frames
    for fn=1:numFrames
        lb(fn, :, :, :) = x(fn, :, :, :) - epsilon;
        ub(fn, :, :, :) = x(fn, :, :, :) + epsilon;
    end

    % Clip the perturbed values to be between 0-1
    lb_min = zeros(numFrames, 120, 160, 3);
    ub_max = ones(numFrames, 120, 160, 3);
    lb_clip = max(lb, lb_min);
    ub_clip = min(ub, ub_max);

    % Create the volume star
    VS = VolumeStar(lb_clip, ub_clip);
end

