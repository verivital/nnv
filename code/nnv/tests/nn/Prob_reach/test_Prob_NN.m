model_name = 'm2nist_dilated_72iou_24layer.mat';
load(model_name)
layers = net.Layers;
newLayers = layers(1:end-2);
Net = SeriesNetwork(newLayers);
obj = matlab2nnv(Net);

load("m2nist_6484_test_images.mat")
image_number = 1;
img = im_data(:,:,image_number);
delta_rgb = 1;
im_lb = img - delta_rgb;
im_ub = img + delta_rgb;
IS = ImageStar(im_lb , im_ub);

reachOptions.train_epochs = 200;
reachOptions.train_lr = 0.0001;
reachOptions.device = "cpu";
reachOptions.train_mode = "cpu";
Reachableset = obj.reachProb_ImageStar(IS, reachOptions);
