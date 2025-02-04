%% Visualize set creation using bias field, l_infinity and adjust contrast

%% Load data

% For this tutorial, choose only one slice
load("data/slice_ISBI_1_1_75.mat");

% Show data
mi_f = min(flair, [], 'all');
ma_f = max(flair,[], 'all');
overlay = labeloverlay(flair,mask,'transparency',0.3);
imshow(overlay,[mi_f, ma_f]);


%% Bias Field

order = 3;
coeff = 0.5;

% img_bf = BiasField(flair, order, coeff);

mi = min(flair, [], 'all');
ma = max(img_bf,[], 'all');

f = figure('Position', get(0, 'Screensize'));
subplot(1,3,1)
imshow(flair, [mi, ma])
colorbar
title("Flair")

subplot(1,3,2)
imshow(img_bf, [mi,ma])
colorbar
title("Flair with Bias Field applied")

im_diff = img_bf-flair;
subplot(1,3,3)
imshow(im_diff, [mi, ma])
colorbar
title("Bias Field difference")



%% Adjust contrast variables

gamma = 0.5;

img_ac = AdjustContrast(flair, gamma); % generates all possible patches to analyze

ma = max(img_ac,[], 'all');

f = figure('Position', get(0, 'Screensize'));
subplot(1,3,1)
imshow(flair, [mi, ma])
colorbar
title("Flair")

subplot(1,3,2)
imshow(img_ac, [mi, ma])
colorbar
title("gammaCorrection(flair,\gamma)")

im_diff = img_ac-flair;
subplot(1,3,3)
imshow(im_diff, [mi, ma])
colorbar
title("UB - LB")



%% L_inf variables

epsilon = 0.2;
nPix = 10; % percetnage of pixels to pertrub

[lb_l, ub_l, ~] = L_inf_transform(flair, wm_mask, epsilon, nPix);

mi = min(lb_l, [], 'all');
ma = max(ub_l,[], 'all');

f = figure('Position', get(0, 'Screensize'));
subplot(1,4,1)
imshow(flair, [mi, ma])
colorbar
title("Flair")

subplot(1,4,2)
imshow(lb_l, [mi,ma])
colorbar
title("LB (Flair - \epsilon)")

subplot(1,4,3)
imshow(ub_l, [mi,ma])
colorbar
title("UB (Flair + \epsilon)")

im_diff = ub_l-lb_l;

subplot(1,4,4)
imshow(im_diff, [min(im_diff, [], 'all'), max(im_diff,[], 'all')])
colorbar
title("UB - LB")