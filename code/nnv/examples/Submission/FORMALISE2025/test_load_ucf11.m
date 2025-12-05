modelName = sprintf("ucf11_c3d_%df_reducedparams_8outchannels.onnx", 16);
% netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "TBCSS", "OutputDataFormats", "BC");
netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "BCSSS", "OutputDataFormats", "BC");
net = matlab2nnv(netonnx);
net.OutputSize = 11;
disp("Finished loading model: " + modelName);

data = readNPY(sprintf("data/UCF11/ucf11_grayscale_%df_verification_data.npy", 16)); % 100 16 112 112 3
labels = readNPY(sprintf("data/UCF11/ucf11_grayscale_%df_verification_labels.npy", 16));
s = data(1,:,:,:,:);
s = squeeze(s);

output = net.evaluate(s);

% TESTED AND OUTPUTS MATCH
%% 

disp(size(s)); % 8 120 160 3
VS = L_inf_attack(s, 1/255, 16);

% Verification settings
verAlg = "relax";

reachOptions = struct;
if verAlg == "relax"
    reachOptions.reachMethod = "relax-star-area";
    reachOptions.relaxFactor = 0.5;
elseif verAlg == "approx"
    reachOptions.reachMethod = "approx-star";
end

t = tic;

% NEED THIS HERE SO MET EXISTS
met = verAlg;

try
    % run verification algorithm
    temp = net.verify_robustness(VS, reachOptions, labels(1)+1);
            
catch ME
    met = ME.message;
    temp = -1;
end

res = temp;
time = toc(t);

fprintf("Result: %d, Time: %f", res, time);
fprintf(met);
fprintf("\n");


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
    lb_min = zeros(numFrames, 112, 112, 1);
    ub_max = ones(numFrames, 112, 112, 1);
    lb_clip = max(lb, lb_min);
    ub_clip = min(ub, ub_max);

    % Create the volume star
    VS = VolumeStar(lb_clip, ub_clip);
end

