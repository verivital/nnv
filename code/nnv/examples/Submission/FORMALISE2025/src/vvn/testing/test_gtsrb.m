%% Load things
% Load data
data = readNPY("../../data/GTSRB/test/gtsrbvideo_16f_test_data_seq.npy");
labels = readNPY("../../data/GTSRB/test/gtsrbvideo_test_labels_seq.npy");

% Preprocessing
% from [B D C H W] to [B D H W C]
reshaped_data = permute(data, [1, 2, 4, 5, 3]);
datacopy = reshaped_data(:,:,:,:,:);

% Experimental variables
numClasses = 43;
% n = 10; % Number of images to evaluate per class
% N = n * numClasses; % Total number of samples to evaluate

% Size of attack
epsilon = [1/255; 2/255; 3/255];
nE = length(epsilon);

% Load the model
modelName = "gtsrb_16f.onnx";
netonnx = importONNXNetwork("../../models/" + modelName, "InputDataFormats", "BCSSS", "OutputDataFormats", "BC");
net = matlab2nnv(netonnx);
net.OutputSize = numClasses;
disp("Finished loading model: " + modelName);

%% Verification settings
reachOptions = struct;
reachOptions.reachMethod = "relax-star-area";
reachOptions.relaxFactor = 0.5;

%% Make predictions on test set to check that we are verifying correct
outputLabels = zeros(length(datacopy));

s = datacopy(1,:,:,:,:);
s = squeeze(s);
l = labels(1) + 1;

output = net.evaluate(s);
[~, P] = max(output);

%%
%%%%%%%%%%%%%%%%
% VERIFICATION %
%%%%%%%%%%%%%%%%

eps = epsilon(1);
fprintf('Starting verification with epsilon %d \n', eps);

% Get the sample
sample = squeeze(datacopy(1,:,:,:,:));

% Perform L_inf attack
VS = L_inf_attack(sample, eps, 16);

%%
t = tic;

% NEED THIS HERE SO MET EXISTS
try
    % run verification algorithm
    fprintf("Verification algorithm starting.")
    temp = net.verify_robustness(VS, reachOptions, labels(1)+1);
    fprintf("Verification algorithm finished.")
            
catch ME
    met = ME.message;
    temp = -1;
    fprintf(met);
end

res = temp;
time = toc(t);

fprintf("\n");
fprintf("Result : %d \n", res);
fprintf("Time: %f", time);

%% Helper Functions
function VS = L_inf_attack(x, epsilon, numFrames)
    lb = squeeze(x);
    ub = squeeze(x);

    % Perturb the frames
    for fn=1:numFrames
        lb(fn, :, :, :) = x(fn, :, :, :) - epsilon;
        ub(fn, :, :, :) = x(fn, :, :, :) + epsilon;
    end

    % Reshape for conversion to VolumeStar
    lb = permute(lb, [2 3 1 4]);
    ub = permute(ub, [2 3 1 4]);

    % Clip the perturbed values to be between 0-1
    lb_min = zeros(32, 32, numFrames, 3);
    ub_max = ones(32, 32, numFrames, 3);
    lb_clip = max(lb, lb_min);
    ub_clip = min(ub, ub_max);

    % Create the volume star
    VS = VolumeStar(lb_clip, ub_clip);
end

