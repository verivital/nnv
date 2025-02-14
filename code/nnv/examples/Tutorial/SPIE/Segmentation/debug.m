% Tutorial example for random noise (L_inf)

% Load one slice from ISBI
load("data/slice_ISBI_1_1_75.mat");

% Load Unet
onnx = importNetworkFromONNX("models/model64.onnx");
% analyzeNetwork(onnx)
net = matlab2nnv(onnx);
windowSize = 64; % from trained model

% Crop background before patches (for all data)
flair = flair(31:158, 31:158);
img_bf = img_bf(31:158, 31:158);
mask = mask(31:158, 31:158);
wm_mask = wm_mask(31:158, 31:158);

% L_inf transform (random noise)
epsilon = 0.004;
nPix = 10; % percetnage of pixels to pertrub
[lb, ub, ~] = L_inf_transform(flair, wm_mask, epsilon, nPix);

% From here, to verify slice, we need to verify 4 patches

% Define verification parameters (reach method)
reachOptions.reachMethod = "relax-star-range-reduceMem";
reachOptions.relaxFactor = 1;
% reachOptions.lp_solver = "gurobi"; % (optional)

img = flair(1:64,1:64);

% Patch 1: Left-top corner
% InputSet1 = ImageStar(lb(1:64,1:64), ub(1:64,1:64));
InputSet1 = ImageStar(img,img);

% Compute reachability
tic;
OutputSet1 = net.reach(InputSet1, reachOptions);
toc;

% Verified output
% verOut1 = verify_output(OutputSet1);

% Predict output
pred1 = onnx.predict(img);

%% Now let's go layer by layer to see where the discrepancy is

nL = length(onnx.Layers);

layerNames = [];
for i = 1:nL
    layerNames = [layerNames; string(onnx.Layers(i).Name)];
end

out = cell(nL, 1);
img = dlarray(img, "SSBC");
for i = 1:nL
    out{i} = forward(onnx, img, 'Outputs', layerNames(i));
end

lbs = cell(nL, 1);
ubs = cell(nL, 1);
for i = 1:nL
    [lbs{i}, ubs{i}] = net.reachSet{i}.estimateRanges();
end

flags = zeros(nL, 1)-1;
be = 1e-4;
% check if every out is inside the range
for i = 1:nL
    flags(i) = all(out{i}+be >= lbs{i},'all') && all(out{i}-be <= ubs{i},'all');
end

