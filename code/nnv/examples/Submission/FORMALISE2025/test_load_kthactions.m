modelName = sprintf("kthactions_%df.onnx", 8);
% netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "TBCSS", "OutputDataFormats", "BC");
netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "BCSSS", "OutputDataFormats", "BC");
net = matlab2nnv(netonnx);
net.OutputSize = 6;
disp("Finished loading model: " + modelName);

data = readNPY(sprintf("data/KTHActions/kthactions_%df.npy", 8)); % 100 8 120 160 3
% reshaped_data = permute(data, [1, 3, 4, 5, 2]); % to match BCSSS %  100 120 160 3 8
reshaped_data = permute(data, [2, 1, 5, 3, 4]); % 8 100 3 120 160
datacopy = reshaped_data(:,:,:,:,:);

% s = datacopy(1,:,:,:,:);
s = data(1,:,:,:,:);
s = squeeze(s);

output = net.evaluate(s);

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

