% Tutorial example for bias field verification

% Load one slice from ISBI
load("data/slice_ISBI_1_1_75.mat");

% Create input set
lb = flair;
ub = img_bf;

windowSize = 64;

% Crop background before patches
flair = flair(35:158, 30:153);

% From here, to verify slice, we need to verify 4 patches

% Left-top corner
