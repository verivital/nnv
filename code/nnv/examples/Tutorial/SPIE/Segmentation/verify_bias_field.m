% Tutorial example for bias field verification

% Load one slice from ISBI
load("data/slice_ISBI_1_1_75.mat");

% Load Unet
onnx = importNetworkFromONNX("models/model64.onnx");
% analyzeNetwork(onnx)
net = matlab2nnv(onnx);
windowSize = 64; % from trained model

% Bias Field
order = 3;
coeff = 0.01;
img_bf = BiasField(flair, order, coeff);

% Crop background before patches (for all data)
flair = flair(31:158, 31:158);
img_bf = img_bf(31:158, 31:158);
mask = mask(31:158, 31:158);

% From here, to verify slice, we need to verify 4 patches

% Define verification parameters (reach method)
reachOptions.reachMethod = "relax-star-area-reduceMem";
reachOptions.relaxFactor = 1;
reachOptions.lp_solver = "gurobi"; % (optional)

% Create bounds of input set
lb = flair;
ub = img_bf;



% Patch 1: Left-top corner
InputSet1 = ImageStar(lb(1:64,1:64), ub(1:64,1:64));

% Compute reachability
tic;
OutputSet1 = net.reach(InputSet1, reachOptions);
toc;

% Verified output
verOut1 = verify_output(OutputSet1);



% Patch 2: Right-top corner
InputSet2 = ImageStar(lb(1:64,65:128), ub(1:64,65:128));

% Compute reachability
tic;
OutputSet2 = net.reach(InputSet2, reachOptions);
toc;

% Verified output
verOut2 = verify_output(OutputSet2);



% Patch 3: Right-bottom corner
InputSet3 = ImageStar(lb(65:128,65:128), ub(65:128,65:128));

% Compute reachability
tic;
OutputSet3 = net.reach(InputSet3, reachOptions);
toc;

% Verified output
verOut3 = verify_output(OutputSet3);



% Patch 4: Lft-bottom corner
InputSet4 = ImageStar(lb(65:128,1:64), ub(65:128,1:64));

% Compute reachability
tic;
OutputSet4 = net.reach(InputSet4, reachOptions);
toc;

% Verified reachability
verOut4 = verify_output(OutputSet4);

% Output reachability
outputSlice = [verOut1 verOut2; verOut3 verOut4];
save("results/biasfield_output.mat","outputSlice")

imshow(outputSlice, [0,2], colormap=hsv(3))
colorbar('XTickLabel', {'Background', 'Lesion', 'Unknown'}, 'XTick',[0,1,2])

% Verified output
verifiedSlice = output_vs_mask(outputSlice, mask);
imshow(verifiedSlice, [-2,2], colormap=hsv(5))
colorbar('XTickLabel', {'False Negative', 'False Positive', 'Background', 'Lesion', 'Unknown'}, 'XTick',[-2,-1,0,1,2])

% Verified lesion
overlay = labeloverlay(flair,verifiedSlice,'transparency',0.3);
imshow(overlay,[-2 2]);

figure;
subplot(1,4,1);
mi_f = min(flair, [], 'all');
ma_f = max(flair,[], 'all');
overlay = labeloverlay(flair,mask,'transparency',0.3);
imshow(overlay,[mi_f, ma_f]);
title("Label mask")

subplot(1,4,2);
overlay = labeloverlay(flair,verifiedSlice==1,'transparency',0.3);
imshow(overlay,[mi_f, ma_f]);
title('Verified lession')

subplot(1,4,3);
overlay = labeloverlay(flair,verifiedSlice==2,'transparency',0.3);
imshow(overlay,[mi_f, ma_f]);
title("Unknown")

subplot(1,4,4);
overlay = labeloverlay(flair,verifiedSlice==-1,'transparency',0.3);
imshow(overlay,[mi_f, ma_f]);
title("False positives")