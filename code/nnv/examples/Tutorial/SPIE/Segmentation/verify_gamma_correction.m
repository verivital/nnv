% Tutorial example for gamma correction verification

% Load one slice from ISBI
load("data/slice_ISBI_1_1_75.mat");

% Load Unet
onnx = importNetworkFromONNX("models/model64.onnx");
% analyzeNetwork(onnx)
net = matlab2nnv(onnx);
windowSize = 64; % from trained model

% Gamma Correction (Adjust contrast)
gamma = 0.95;
img_ac = AdjustContrast(flair, gamma); % generates all possible patches to analyze

% Crop background before patches (for all data)
flair = flair(31:158, 31:158);
img_ac = img_ac(31:158, 31:158);
mask = mask(31:158, 31:158);

% From here, to verify slice, we need to verify 4 patches

% Define verification parameters (reach method)
reachOptions.reachMethod = "relax-star-area";
reachOptions.relaxFactor = 1;
% reachOptions.lp_solver = "gurobi"; % (optional)

% Create bounds of input set
lb = flair;
ub = img_ac;

% Patch 1: Left-top corner
InputSet1 = ImageStar(lb(1:64,1:64), ub(1:64,1:64));

% Compute reachability
tic;
OutputSet1 = net.reach(InputSet1, reachOptions);
toc;

% Verified output
verOut1 = verify_output(OutputSet1);

pred1 = onnx.predict(flair(1:64,1:64));


% Patch 2: Right-top corner
InputSet2 = ImageStar(lb(1:64,65:128), ub(1:64,65:128));

% Compute reachability
tic;
OutputSet2 = net.reach(InputSet2, reachOptions);
toc;

% Verified output
verOut2 = verify_output(OutputSet2);

pred2 = onnx.predict(flair(1:64,65:128));


% Patch 3: Right-bottom corner
InputSet3 = ImageStar(lb(65:128,65:128), ub(65:128,65:128));

% Compute reachability
tic;
OutputSet3 = net.reach(InputSet3, reachOptions);
toc;

% Verified output
verOut3 = verify_output(OutputSet3);

pred3 = onnx.predict(flair(65:128,65:128));

% Patch 4: Lft-bottom corner
InputSet4 = ImageStar(lb(65:128,1:64), ub(65:128,1:64));

% Compute reachability
tic;
OutputSet4 = net.reach(InputSet4, reachOptions);
toc;

% Verified output
verOut4 = verify_output(OutputSet4);

pred4 = onnx.predict(flair(65:128,1:64));

% Output reachability
outputSlice = [verOut1 verOut2; verOut3 verOut4];
pred = [pred1 pred2; pred3 pred4];
pred = pred > 0;
save("results/gama_output.mat","outputSlice")


%%
custommap = [0.2 0.1 0.5
    0.8 0.7 0.3
    0.8 0 0.2
    0.1 0.5 0.8
    0.2 0.7 0.6
    ];

% figure;
% imshow(pred, [0,1], colormap=custommap(1:2,:))
% colorbar('XTickLabel', {'Background', 'Lesion'}, 'XTick',[0,1])

pred_label = pred_vs_mask(pred,mask);
figure;
imshow(pred_label, [0,4], colormap=custommap([1:2,4:5],:))
colorbar('XTickLabel', {'Background', 'Lesion', 'False Positive', 'False Negative'}, 'XTick',[0,1,3,4])

% figure;
% imshow(outputSlice, [0,2], colormap=custommap(1:3,:))
% colorbar('XTickLabel', {'Background', 'Lesion', 'Unknown'}, 'XTick',[0,1,2])

% Verified output
verifiedSlice = output_vs_mask(outputSlice, mask, pred);
figure;
imshow(verifiedSlice, [0,4], colormap=custommap)
colorbar('XTickLabel', {'Background', 'Lesion', 'Unknown', 'False Positive', 'False Negative'}, 'XTick',[0,1,2,3,4])

% Verified lesion
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
overlay = labeloverlay(flair,verifiedSlice==0,'transparency',0.3);
imshow(overlay,[mi_f, ma_f]);
title("Verified Background")

