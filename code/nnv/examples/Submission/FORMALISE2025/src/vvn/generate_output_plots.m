%% Load things
% Load data
data = readNPY("data/STMNIST/stmnistvideo_64f_test_data_seq.npy");
labels = readNPY("data/STMNIST/stmnistvideo_64f_test_labels_seq.npy");

% Preprocessing
% from [B D C H W] to [B D H W C]
reshaped_data = permute(data, [1, 3, 4, 5, 2]);
datacopy = reshaped_data(:,:,:,:,:);

% Experimental variables
numClasses = 10;

% Size of attack
epsilon = [1/255; 2/255; 3/255];
nE = length(epsilon);

% Load the model
modelName = "stmnist_64f.onnx";
netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "BCSSS", "OutputDataFormats", "BC");
net = matlab2nnv(netonnx);
net.OutputSize = numClasses;
disp("Finished loading model: " + modelName);

%% Verification settings
reachOptions = struct;
reachOptions.reachMethod = "relax-star-area";
reachOptions.relaxFactor = 0.5;

%% Make predictions on test set to check that we are verifying correct
outputLabels = zeros(length(datacopy));
index = 3; % 311 for a really long example

s = datacopy(index,:,:,:,:);
s = squeeze(s);
l = labels(index) + 1;

fprintf('Label: %d \n', l);

%% ROBUST SAMPLE
eps = 1/255;

fprintf('Starting verification with epsilon %d \n', eps);

% Perform L_inf attack
[VS, lb_clip, ub_clip] = L_inf_attack(s, eps, 64);

output = net.evaluate(s);
[~, P] = max(output);

LB_output = net.evaluate(lb_clip);
[~, LBPred] = max(LB_output);

UB_output = net.evaluate(ub_clip);
[~, UBPred] = max(UB_output);

%%
t = tic;

% NEED THIS HERE SO MET EXISTS
try
    % run verification algorithm
    fprintf("Verification algorithm starting.\n")
    temp = net.verify_robustness(VS, reachOptions, l);
    fprintf("Verification algorithm finished.\n")
            
catch ME
    met = ME.message;
    temp = -1;
    fprintf(met);
end

p = profile('info');

res = temp;
time = toc(t);

fprintf("\n");
fprintf("Result : %d \n", res);
fprintf("Time: %f\n", time);

fprintf("Computing reach sets...\n")
% Get the reachable sets 
R = net.reachSet{end};

fprintf("Done computing reach sets! \n")


fprintf("Get the ranges for ub/lb \n")
% Get ranges for each output index
[lb_out, ub_out] = R.getRanges;
lb_out = squeeze(lb_out);
ub_out = squeeze(ub_out);

fprintf("Now to plotting! \n");

% Get middle point for each output and range sizes
mid_range = (lb_out + ub_out)/2;
range_size = ub_out - mid_range;

% Label for x-axis
x = [0 1 2 3 4 5 6 7 8 9];

% Visualize set ranges and evaluation points
f = figure;
errorbar(x, mid_range, range_size, '.', 'Color', 'r', 'LineWidth', 2);
hold on;
xlim([-0.5, 9.5]);
scatter(x, output, 30, 'x', 'MarkerEdgeColor', 'r');
title('Reachable Outputs');
xlabel('Label');
ylabel('Reachable Output Range on the Input Set');
% Save the figure
saveas(f, "/tmp/reach_stmnist_plot.png");

%% NON-ROBUST SAMPLE
eps = 3/255;

fprintf('Starting verification with epsilon %d \n', eps);

% Perform L_inf attack
[VS, lb_clip, ub_clip] = L_inf_attack(s, eps, 64);

output = net.evaluate(s);
[~, P] = max(output);

LB_output = net.evaluate(lb_clip);
[~, LBPred] = max(LB_output);

UB_output = net.evaluate(ub_clip);
[~, UBPred] = max(UB_output);

%%
t = tic;

% NEED THIS HERE SO MET EXISTS
try
    % run verification algorithm
    fprintf("Verification algorithm starting.\n")
    temp = net.verify_robustness(VS, reachOptions, l);
    fprintf("Verification algorithm finished.\n")
            
catch ME
    met = ME.message;
    temp = -1;
    fprintf(met);
end

p = profile('info');

res = temp;
time = toc(t);

fprintf("\n");
fprintf("Result : %d \n", res);
fprintf("Time: %f\n", time);

fprintf("Computing reach sets...\n")
% Get the reachable sets 
R = net.reachSet{end};

fprintf("Done computing reach sets! \n")


fprintf("Get the ranges for ub/lb \n")
% Get ranges for each output index
[lb_out, ub_out] = R.getRanges;
lb_out = squeeze(lb_out);
ub_out = squeeze(ub_out);

fprintf("Now to plotting! \n");

% Get middle point for each output and range sizes
mid_range = (lb_out + ub_out)/2;
range_size = ub_out - mid_range;

% Label for x-axis
x = [0 1 2 3 4 5 6 7 8 9];

% Visualize set ranges and evaluation points
f = figure;
errorbar(x, mid_range, range_size, '.', 'Color', 'r', 'LineWidth', 2);
hold on;
xlim([-0.5, 9.5]);
scatter(x, output, 30, 'x', 'MarkerEdgeColor', 'r');
title('Reachable Outputs');
xlabel('Label');
ylabel('Reachable Output Range on the Input Set');
% Save the figure
saveas(f, "/tmp/reach_bad_stmnist_plot.png");


%% Helper Functions
function [VS, lb_clip, ub_clip] = L_inf_attack(x, epsilon, numFrames)
    lb = squeeze(x);
    ub = squeeze(x);

    % Perturb the frames
    for fn=1:numFrames
        lb(fn, :, :, :) = x(fn, :, :, :) - epsilon;
        ub(fn, :, :, :) = x(fn, :, :, :) + epsilon;
    end

    % Clip the perturbed values to be between 0-1
    lb_min = zeros(numFrames, 10, 10, 2);
    ub_max = ones(numFrames, 10, 10, 2);
    lb_clip = max(lb, lb_min);
    ub_clip = min(ub, ub_max);

    % Create the volume star
    VS = VolumeStar(lb_clip, ub_clip);
end

