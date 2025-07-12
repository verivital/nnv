function Reachableset = test_Prob_utils()

model_name = 'm2nist_dilated_72iou_24layer.mat';
dir0 = pwd;
parentDir = fileparts(dir0);
parentDir = fileparts(parentDir);
nnv_dir = fileparts(parentDir);
addpath(genpath(nnv_dir))

load(model_name)
layers = net.Layers;
newLayers = layers(1:end-2);
Net = SeriesNetwork(newLayers);

load("m2nist_6484_test_images.mat")
image_number = 1;
img = im_data(:,:,image_number);

delta_rgb = 1;
height = 64;
width = 84;
n_channel = 1;
train_epochs = 200;
train_lr = 0.0001;

im_lb = img - delta_rgb;
im_ub = img + delta_rgb;

ct = 0;
indices = zeros(height*width*n_channel , 2);
for i=1:height
    for j=1:width
        ct =ct+1;
        indices(ct , :) = [i,j];
    end
end


reachOptions.train_epochs = train_epochs;
reachOptions.train_lr = train_lr;
reachOptions.dims = [-1 -1];
reachOptions.coverage = 0.99;
reachOptions.confidence = 0.9;
reachOptions.train_mode = 'Linear';
reachOptions.surrogate_dim = [];
reachOptions.threshold_normal = 1e-5;
reachOptions.indices = indices;



IS = ImageStar(im_lb , im_ub);

Reachableset  = Prob_reach(Net, IS, reachOptions);


end